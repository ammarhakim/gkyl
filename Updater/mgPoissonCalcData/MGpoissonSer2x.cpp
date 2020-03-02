#include <MGpoissonModDecl.h> 
 
void MGpoissonProlong2xSer_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF1 = fldF[0];
  double *fldF2 = fldF[1];
  double *fldF3 = fldF[2];
  double *fldF4 = fldF[3];

  fldF1[0] = 0.75*fldC[3]-0.8660254037844386*fldC[2]-0.8660254037844386*fldC[1]+fldC[0]; 
  fldF1[1] = 0.5*fldC[1]-0.4330127018922193*fldC[3]; 
  fldF1[2] = 0.5*fldC[2]-0.4330127018922193*fldC[3]; 
  fldF1[3] = 0.25*fldC[3]; 

  fldF2[0] = (-0.75*fldC[3])-0.8660254037844386*fldC[2]+0.8660254037844386*fldC[1]+fldC[0]; 
  fldF2[1] = 0.5*fldC[1]-0.4330127018922193*fldC[3]; 
  fldF2[2] = 0.4330127018922193*fldC[3]+0.5*fldC[2]; 
  fldF2[3] = 0.25*fldC[3]; 

  fldF3[0] = (-0.75*fldC[3])+0.8660254037844386*fldC[2]-0.8660254037844386*fldC[1]+fldC[0]; 
  fldF3[1] = 0.4330127018922193*fldC[3]+0.5*fldC[1]; 
  fldF3[2] = 0.5*fldC[2]-0.4330127018922193*fldC[3]; 
  fldF3[3] = 0.25*fldC[3]; 

  fldF4[0] = 0.75*fldC[3]+0.8660254037844386*fldC[2]+0.8660254037844386*fldC[1]+fldC[0]; 
  fldF4[1] = 0.4330127018922193*fldC[3]+0.5*fldC[1]; 
  fldF4[2] = 0.4330127018922193*fldC[3]+0.5*fldC[2]; 
  fldF4[3] = 0.25*fldC[3]; 

}

void MGpoissonRestrict2xSer_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in stencils pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[0];
  double *fldF2 = fldF[1];
  double *fldF3 = fldF[2];
  double *fldF4 = fldF[3];

  fldC[0] = 0.25*fldF4[0]+0.25*fldF3[0]+0.25*fldF2[0]+0.25*fldF1[0]; 
  fldC[1] = 0.125*fldF4[1]+0.125*fldF3[1]+0.125*fldF2[1]+0.125*fldF1[1]+0.2165063509461096*fldF4[0]-0.2165063509461096*fldF3[0]+0.2165063509461096*fldF2[0]-0.2165063509461096*fldF1[0]; 
  fldC[2] = 0.125*fldF4[2]+0.125*fldF3[2]+0.125*fldF2[2]+0.125*fldF1[2]+0.2165063509461096*fldF4[0]+0.2165063509461096*fldF3[0]-0.2165063509461096*fldF2[0]-0.2165063509461096*fldF1[0]; 
  fldC[3] = 0.0625*fldF4[3]+0.0625*fldF3[3]+0.0625*fldF2[3]+0.0625*fldF1[3]+0.1082531754730548*fldF4[2]-0.1082531754730548*fldF3[2]+0.1082531754730548*fldF2[2]-0.1082531754730548*fldF1[2]+0.1082531754730548*fldF4[1]+0.1082531754730548*fldF3[1]-0.1082531754730548*fldF2[1]-0.1082531754730548*fldF1[1]+0.1875*fldF4[0]-0.1875*fldF3[0]-0.1875*fldF2[0]+0.1875*fldF1[0]; 

}

void MGpoissonGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0])*rdx2SqVol[0])/(18.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[1] = (16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])/(18.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 
  phiC[2] = (16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]+8.660254037844386*rdx2SqVol[0]*phiLx[3]-7.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]-7.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1])/(46.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[3] = (16.0*rho[3]*volFac-7.0*rdx2SqVol[1]*phiUy[3]-7.0*rdx2SqVol[0]*phiUx[3]-7.0*rdx2SqVol[1]*phiLy[3]-7.0*rdx2SqVol[0]*phiLx[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1])/(46.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[0]*rho[1]+288.0*rho[0]*rdx2SqVol[1]+2096.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(162.0*phiUy[0]+162.0*phiLy[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]-31.17691453623978*rdx2SqVol[0]*phiUx[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]+(1179.0*phiUy[0]+54.0*phiUx[0]+1179.0*phiLy[0]+1296.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiUx[1]+(1416.0*phiUx[0]+4224.0*bcVals[0])*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]+471.1178196587346*rdx2SqVol[0]*rho[0])*volFac+((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]+162.0*phiLy[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]-342.0*rdx2SqVol[0]*phiUx[1]+351.0*rdx2SqVol[0]*phiLy[1]+(265.0037735580381*phiUy[0]+342.9460598986376*phiUx[0]+265.0037735580381*phiLy[0]-1745.907214029428*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiUx[1]+(831.384387633061*phiUx[0]-1662.768775266122*bcVals[0])*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = ((859.0972005541631*rdx2SqVol[0]*rho[3]+(736.0*rdx2SqVol[1]+2096.0*rdx2SqVol[0])*rho[2])*volFac-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(138.0*rdx2SqVol[0]*rdx2SqVol[1]+1416.0*rdx2SqVolSq[0])*phiUx[2]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(398.3716857408418*phiUy[0]-398.3716857408418*phiLy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiUy[0]-1134.493278957615*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[0]*rho[2])*volFac+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiUx[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+(876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1]+831.384387633061*rdx2SqVolSq[0])*phiUx[2]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiUy[0]-255.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[1]-864.0*rho[0]*rdx2SqVol[1]-2256.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-486.0*phiUy[0])-486.0*phiLy[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+467.6537180435967*rdx2SqVol[0]*phiUx[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+((-1269.0*phiUy[0])-486.0*phiUx[0]-1269.0*phiLy[0]+864.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiUx[1]+(2816.0*bcVals[0]-984.0*phiUx[0])*rdx2SqVolSq[0]))/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]-526.5434455009387*rdx2SqVol[0]*rho[0])*volFac+((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]+486.0*phiLy[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]-522.0*rdx2SqVol[0]*phiUx[1]+243.0*rdx2SqVol[0]*phiLy[1]+((-296.1806880942779*phiUy[0])+592.3613761885558*phiUx[0]-296.1806880942779*phiLy[0]+1163.938142686285*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiUx[1]+1108.512516844081*bcVals[0]*rdx2SqVolSq[0])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[3]+((-2208.0*rdx2SqVol[1])-2256.0*rdx2SqVol[0])*rho[2])*volFac-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-1242.0*rdx2SqVol[0]*rdx2SqVol[1])-984.0*rdx2SqVolSq[0])*phiUx[2]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1195.115057222525*phiLy[0]-1195.115057222525*phiUy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiLy[0]-1221.095819336058*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1]))/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[0]*rho[2])*volFac+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUx[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiLy[0]-285.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[1]-288.0*rho[0]*rdx2SqVol[1]-2096.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-162.0*phiUy[0])-162.0*phiLy[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]-31.17691453623978*rdx2SqVol[0]*phiLx[1]-1296.0*rdx2SqVol[0]*bcVals[1]+((-1179.0*phiUy[0])-1179.0*phiLy[0]-54.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiLx[1]-4224.0*rdx2SqVolSq[0]*bcVals[1]-1416.0*phiLx[0]*rdx2SqVolSq[0]))/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]-471.1178196587346*rdx2SqVol[0]*rho[0])*volFac+((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]+162.0*phiLy[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]+351.0*rdx2SqVol[0]*phiLy[1]-342.0*rdx2SqVol[0]*phiLx[1]+1745.907214029428*rdx2SqVol[0]*bcVals[1]+((-265.0037735580381*phiUy[0])-265.0037735580381*phiLy[0]-342.9460598986376*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiLx[1]+1662.768775266122*rdx2SqVolSq[0]*bcVals[1]-831.384387633061*phiLx[0]*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[3]+((-736.0*rdx2SqVol[1])-2096.0*rdx2SqVol[0])*rho[2])*volFac-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiLx[3]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-138.0*rdx2SqVol[0]*rdx2SqVol[1])-1416.0*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiLy[0]-398.3716857408418*phiUy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiLy[0]-1134.493278957615*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1]))/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[0]*rho[2])*volFac+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiLx[3]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+((-876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1])-831.384387633061*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiLy[0]-255.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[0]*rho[1]+864.0*rho[0]*rdx2SqVol[1]+2256.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(486.0*phiUy[0]+486.0*phiLy[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+467.6537180435967*rdx2SqVol[0]*phiLx[1]+864.0*rdx2SqVol[0]*bcVals[1]+(1269.0*phiUy[0]+1269.0*phiLy[0]+486.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiLx[1]+2816.0*rdx2SqVolSq[0]*bcVals[1]+984.0*phiLx[0]*rdx2SqVolSq[0])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]+526.5434455009387*rdx2SqVol[0]*rho[0])*volFac+((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]+486.0*phiLy[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]+243.0*rdx2SqVol[0]*phiLy[1]-522.0*rdx2SqVol[0]*phiLx[1]+1163.938142686285*rdx2SqVol[0]*bcVals[1]+(296.1806880942779*phiUy[0]+296.1806880942779*phiLy[0]-592.3613761885558*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiLx[1]+1108.512516844081*rdx2SqVolSq[0]*bcVals[1])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = ((415.6921938165305*rdx2SqVol[0]*rho[3]+(2208.0*rdx2SqVol[1]+2256.0*rdx2SqVol[0])*rho[2])*volFac-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiLx[3]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1242.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiLx[2]+(1195.115057222525*phiUy[0]-1195.115057222525*phiLy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiUy[0]-1221.095819336058*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[0]*rho[2])*volFac+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLx[3]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]-1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiUy[0]-285.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[1]*rho[2]+2096.0*rho[0]*rdx2SqVol[1]+288.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(4224.0*rdx2SqVolSq[1]+1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+1416.0*phiUy[0]*rdx2SqVolSq[1]+((-1134.493278957615*rdx2SqVol[0]*phiUx[1])+1134.493278957615*rdx2SqVol[0]*phiLx[1]+(54.0*phiUy[0]+1179.0*phiUx[0]+1179.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-155.8845726811989*rdx2SqVolSq[0]*phiUx[1]+155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+(162.0*phiUx[0]+162.0*phiLx[0])*rdx2SqVolSq[0])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = ((859.0972005541631*rdx2SqVol[1]*rho[3]+(2096.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[1])*volFac+((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+1416.0*phiUy[1]*rdx2SqVolSq[1]+(138.0*rdx2SqVol[0]*phiUy[1]-917.0*rdx2SqVol[0]*phiUx[1]-917.0*rdx2SqVol[0]*phiLx[1]+(1134.493278957615*phiUx[0]-1134.493278957615*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-322.0*rdx2SqVolSq[0]*phiUx[1]-322.0*rdx2SqVolSq[0]*phiLx[1]+(398.3716857408418*phiUx[0]-398.3716857408418*phiLx[0])*rdx2SqVolSq[0])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]+471.1178196587346*rho[0]*rdx2SqVol[1])*volFac+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]+((-1662.768775266122*rdx2SqVolSq[1])-1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+831.384387633061*phiUy[0]*rdx2SqVolSq[1]+((-255.0*rdx2SqVol[0]*phiUx[1])+255.0*rdx2SqVol[0]*phiLx[1]+(342.9460598986376*phiUy[0]+265.0037735580381*phiUx[0]+265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[1]*rho[1])*volFac+((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]+831.384387633061*phiUy[1]*rdx2SqVolSq[1]+(876.4177086298519*rdx2SqVol[0]*phiUy[1]-206.1140461006964*rdx2SqVol[0]*phiUx[1]-206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiUx[0]-255.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[2]-2256.0*rho[0]*rdx2SqVol[1]-864.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]-984.0*phiUy[0]*rdx2SqVolSq[1]+(1221.095819336058*rdx2SqVol[0]*phiUx[1]-1221.095819336058*rdx2SqVol[0]*phiLx[1]+((-486.0*phiUy[0])-1269.0*phiUx[0]-1269.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0]*phiUx[1]-467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+((-486.0*phiUx[0])-486.0*phiLx[0])*rdx2SqVolSq[0]))/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[3]+((-2256.0*rdx2SqVol[1])-2208.0*rdx2SqVol[0])*rho[1])*volFac+(969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]-984.0*phiUy[1]*rdx2SqVolSq[1]+((-1242.0*rdx2SqVol[0]*phiUy[1])+987.0*rdx2SqVol[0]*phiUx[1]+987.0*rdx2SqVol[0]*phiLx[1]+(1221.095819336058*phiLx[0]-1221.095819336058*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+966.0*rdx2SqVolSq[0]*phiUx[1]+966.0*rdx2SqVolSq[0]*phiLx[1]+(1195.115057222525*phiLx[0]-1195.115057222525*phiUx[0])*rdx2SqVolSq[0]))/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]-526.5434455009387*rho[0]*rdx2SqVol[1])*volFac+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(285.0*rdx2SqVol[0]*phiUx[1]-285.0*rdx2SqVol[0]*phiLx[1]+(592.3613761885558*phiUy[0]-296.1806880942779*phiUx[0]-296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[1]*rho[1])*volFac+(24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+(1513.812405815199*rdx2SqVol[0]*phiUy[1]+230.3627574066607*rdx2SqVol[0]*phiUx[1]+230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiLx[0]-285.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[2]-2096.0*rho[0]*rdx2SqVol[1]-288.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-4224.0*rdx2SqVolSq[1])-1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]-1416.0*phiLy[0]*rdx2SqVolSq[1]+(1134.493278957615*rdx2SqVol[0]*phiUx[1]-1134.493278957615*rdx2SqVol[0]*phiLx[1]+((-1179.0*phiUx[0])-54.0*phiLy[0]-1179.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0]*phiUx[1]-155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+((-162.0*phiUx[0])-162.0*phiLx[0])*rdx2SqVolSq[0]))/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[3]+((-2096.0*rdx2SqVol[1])-736.0*rdx2SqVol[0])*rho[1])*volFac-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]-1416.0*phiLy[1]*rdx2SqVolSq[1]+(917.0*rdx2SqVol[0]*phiUx[1]-138.0*rdx2SqVol[0]*phiLy[1]+917.0*rdx2SqVol[0]*phiLx[1]+(1134.493278957615*phiLx[0]-1134.493278957615*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+322.0*rdx2SqVolSq[0]*phiUx[1]+322.0*rdx2SqVolSq[0]*phiLx[1]+(398.3716857408418*phiLx[0]-398.3716857408418*phiUx[0])*rdx2SqVolSq[0]))/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]-471.1178196587346*rho[0]*rdx2SqVol[1])*volFac+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+(1662.768775266122*rdx2SqVolSq[1]+1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[0]*rdx2SqVolSq[1]+(255.0*rdx2SqVol[0]*phiUx[1]-255.0*rdx2SqVol[0]*phiLx[1]+((-265.0037735580381*phiUx[0])-342.9460598986376*phiLy[0]-265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[1]*rho[1])*volFac+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[1]*rdx2SqVolSq[1]+(206.1140461006964*rdx2SqVol[0]*phiUx[1]-876.4177086298519*rdx2SqVol[0]*phiLy[1]+206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiLx[0]-255.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[1]*rho[2]+2256.0*rho[0]*rdx2SqVol[1]+864.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+984.0*phiLy[0]*rdx2SqVolSq[1]+((-1221.095819336058*rdx2SqVol[0]*phiUx[1])+1221.095819336058*rdx2SqVol[0]*phiLx[1]+(1269.0*phiUx[0]+486.0*phiLy[0]+1269.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-467.6537180435967*rdx2SqVolSq[0]*phiUx[1]+467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+(486.0*phiUx[0]+486.0*phiLx[0])*rdx2SqVolSq[0])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = ((415.6921938165305*rdx2SqVol[1]*rho[3]+(2256.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[1])*volFac-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+(969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+984.0*phiLy[1]*rdx2SqVolSq[1]+((-987.0*rdx2SqVol[0]*phiUx[1])+1242.0*rdx2SqVol[0]*phiLy[1]-987.0*rdx2SqVol[0]*phiLx[1]+(1221.095819336058*phiUx[0]-1221.095819336058*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-966.0*rdx2SqVolSq[0]*phiUx[1]-966.0*rdx2SqVolSq[0]*phiLx[1]+(1195.115057222525*phiUx[0]-1195.115057222525*phiLx[0])*rdx2SqVolSq[0])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]+526.5434455009387*rho[0]*rdx2SqVol[1])*volFac+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+((-285.0*rdx2SqVol[0]*phiUx[1])+285.0*rdx2SqVol[0]*phiLx[1]+(296.1806880942779*phiUx[0]-592.3613761885558*phiLy[0]+296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[1]*rho[1])*volFac+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+(24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+((-230.3627574066607*rdx2SqVol[0]*phiUx[1])-1513.812405815199*rdx2SqVol[0]*phiLy[1]-230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiUx[0]-285.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+624456.0*phiUy[0]*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiUx[1]+(3368931.0*phiUy[0]+173313.0*phiUx[0]+4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiUx[1]+(173313.0*phiUy[0]+3368931.0*phiUx[0]+1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+(624456.0*phiUx[0]+1862784.0*bcVals[0])*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+624456.0*phiUy[1]*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiUx[1]+(847040.3948826761*phiUy[0]+1100685.379244677*phiUx[0]-5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiUx[1]+(38955.5547130316*phiUy[0]+2275410.734360502*phiUx[0]-5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiUx[1]+(366640.5149461798*phiUx[0]-733281.0298923596*bcVals[0])*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiUx[1]+(2275410.734360502*phiUy[0]+38955.5547130316*phiUx[0]+934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiUx[1]+(1100685.379244677*phiUy[0]+847040.3948826761*phiUx[0]+1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiUx[1]+(90933.0*phiUy[0]+82467.0*phiUx[0]-419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiUx[1]+(82467.0*phiUy[0]+90933.0*phiUx[0]+73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(17043.37994647775*rdx2SqVolCu[1]+381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]+973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1])-2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]-1136585.596333157*rdx2SqVolCu[0])*rho[1]-92496.0*rho[0]*rdx2SqVolCu[1]-2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0]*rho[0])*volFac+(307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(39767.88654178142*rdx2SqVolR4[1]+930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]+2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1]+257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]+747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(115456.0*rdx2SqVolR4[1]+2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]+7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-40344.0*phiUy[0]*rdx2SqVolR4[1]+((-310402.557275226*rdx2SqVol[0]*phiUy[1])+10012.98571855568*rdx2SqVol[0]*phiUx[1]+((-945501.0*phiUy[0])-17343.0*phiUx[0]-416232.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1122732.653974222*rdx2SqVolSq[0]*phiUy[1])+1113691.348758712*rdx2SqVolSq[0]*phiUx[1]+((-2972058.0*phiUy[0])-1286154.0*phiUx[0]-5155056.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-639329.3979374008*rdx2SqVolCu[0]*phiUy[1])+3757176.73613406*rdx2SqVolCu[0]*phiUx[1]+((-1559817.0*phiUy[0])-4278411.0*phiUx[0]-1.3513464e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1649882.317257809*rdx2SqVolR4[0]*phiUx[1]+((-1873368.0*phiUx[0])-5588352.0*bcVals[0])*rdx2SqVolR4[0]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(85680.0*rdx2SqVol[0]*rdx2SqVolSq[1]+260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-92496.0*rdx2SqVolCu[1])-873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]-825552.0*rdx2SqVolCu[0])*rho[1]-407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0]*rho[0])*volFac+(39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(168555.0*rdx2SqVol[0]*rdx2SqVolCu[1]+606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(20295.0*rdx2SqVol[0]*rdx2SqVolCu[1]+151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1]+1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]+623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-40344.0*phiUy[1]*rdx2SqVolR4[1]+((-413649.0*rdx2SqVol[0]*phiUy[1])+109839.0*rdx2SqVol[0]*phiUx[1]+((-170220.7572154465*phiUy[0])-110142.8429041125*phiUx[0]+560727.200239118*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1048338.0*rdx2SqVolSq[0]*phiUy[1])+1081578.0*rdx2SqVolSq[0]*phiUx[1]+((-615692.1005665087*phiUy[0])-1116705.117163882*phiUx[0]+3464794.435460782*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-464373.0*rdx2SqVolCu[0]*phiUy[1])+2599263.0*rdx2SqVolCu[0]*phiUx[1]+((-350599.9924172843*phiUy[0])-2717894.290068507*phiUx[0]+6136988.564971584*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1047816.0*rdx2SqVolR4[0]*phiUx[1]+(2199843.089677078*bcVals[0]-1099921.544838539*phiUx[0])*rdx2SqVolR4[0]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = (((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+(17712.0*rdx2SqVolCu[1]+472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0])*rho[2]+((-197904.0*rdx2SqVol[0]*rdx2SqVolSq[1])-600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-21588.28126553848*rho[0]*rdx2SqVolCu[1]-482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiUx[3]+(984.0*rdx2SqVolR4[1]-24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]-659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(3321.0*rdx2SqVol[0]*rdx2SqVolCu[1]+156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]+1873368.0*rdx2SqVolR4[0])*phiUx[2]+(45449.01319060734*rdx2SqVolR4[1]+1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]+4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]+3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-72447.0*rdx2SqVol[0]*phiUy[1])+2337.0*rdx2SqVol[0]*phiUx[1]+(52621.43558475006*phiUy[0]-4047.802737288466*phiUx[0]-97147.26569492316*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(287280.0*rdx2SqVolSq[0]*phiUx[1]+(812719.8081306986*phiUy[0]-326193.6644878315*phiUx[0]-973052.2872857346*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(779247.0*rdx2SqVolCu[0]*phiUy[1]+846963.0*rdx2SqVolCu[0]*phiUx[1]+(1901183.83687717*phiUy[0]-946692.2060453439*phiUx[0]-1908983.261663653*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+(31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1]+270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0])*rho[2]+((-21588.28126553848*rdx2SqVolCu[1])-143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]-204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiUx[3]+(39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1]+15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]-376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1]+453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]+1099921.544838539*rdx2SqVolR4[0])*phiUx[2]+(150416.0*rdx2SqVol[0]*rdx2SqVolCu[1]+693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(176754.0528615963*rdx2SqVol[0]*phiUy[1]+25636.08400282695*rdx2SqVol[0]*phiUx[1]+((-39729.0*phiUy[0])-25707.0*phiUx[0]+130872.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(812719.8081306986*rdx2SqVolSq[0]*phiUy[1]+182447.3038660752*rdx2SqVolSq[0]*phiUx[1]+(383040.0*bcVals[0]-191520.0*phiUx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(566001.2949481651*rdx2SqVolCu[0]*phiUy[1]+278113.666120527*rdx2SqVolCu[0]*phiUx[1]+(427329.0*phiUy[0]-304893.0*phiUx[0]-244872.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-1136585.596333157*rdx2SqVolCu[1])-2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]-742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1]+381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[1]-2773008.0*rho[0]*rdx2SqVolCu[1]-6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0]*rho[0])*volFac+((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(1649882.317257809*rdx2SqVolR4[1]+3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]+1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]+10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1])-1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]-310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+((-5588352.0*rdx2SqVolR4[1])-1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1873368.0*phiUy[0]*rdx2SqVolR4[1]+(747388.5837200081*rdx2SqVol[0]*phiUy[1]+1500934.608060923*rdx2SqVol[0]*phiUx[1]+((-4278411.0*phiUy[0])-1559817.0*phiUx[0]+2773008.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(257521.3140693406*rdx2SqVolSq[0]*phiUy[1]+2913967.637637727*rdx2SqVolSq[0]*phiUx[1]+((-1286154.0*phiUy[0])-2972058.0*phiUx[0]+7782592.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3195.633739964578*rdx2SqVolCu[0]*phiUy[1]+930985.9693223094*rdx2SqVolCu[0]*phiUx[1]+((-17343.0*phiUy[0])-945501.0*phiUx[0]+2659024.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0]*phiUx[1]+(115456.0*bcVals[0]-40344.0*phiUx[0])*rdx2SqVolR4[0]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = (((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-600780.0*rdx2SqVol[0]*rdx2SqVolSq[1])-197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(2773008.0*rdx2SqVolCu[1]+2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]+472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[1]-1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0]*rho[0])*volFac+((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(846963.0*rdx2SqVol[0]*rdx2SqVolCu[1]+287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(779247.0*rdx2SqVol[0]*rdx2SqVolCu[1]-72447.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+((-1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1])-973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]-97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1873368.0*phiUy[1]*rdx2SqVolR4[1]+(998973.0*rdx2SqVol[0]*phiUy[1]-1675359.0*rdx2SqVol[0]*phiUx[1]+((-946692.2060453439*phiUy[0])+1901183.83687717*phiUx[0]+3735659.468951633*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(156186.0*rdx2SqVolSq[0]*phiUy[1]-659958.0*rdx2SqVolSq[0]*phiUx[1]+((-326193.6644878315*phiUy[0])+812719.8081306986*phiUx[0]+4193890.830602055*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3321.0*rdx2SqVolCu[0]*phiUy[1]-24363.0*rdx2SqVolCu[0]*phiUx[1]+((-4047.802737288466*phiUy[0])+52621.43558475006*phiUx[0]+1119902.482954654*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+984.0*rdx2SqVolR4[0]*phiUx[1]+45449.01319060734*bcVals[0]*rdx2SqVolR4[0])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+((-825552.0*rdx2SqVolCu[1])-2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]-873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0])*rho[2]+(260100.0*rdx2SqVol[0]*rdx2SqVolSq[1]+85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-623288.8754085057*rho[0]*rdx2SqVolCu[1]-1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiUx[3]+(1047816.0*rdx2SqVolR4[1]+2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-464373.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiUx[2]+(2199843.089677078*rdx2SqVolR4[1]+6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]+3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]+560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1099921.544838539*phiUy[0]*rdx2SqVolR4[1]+(240705.0*rdx2SqVol[0]*phiUy[1]+337365.0*rdx2SqVol[0]*phiUx[1]+((-2717894.290068507*phiUy[0])-350599.9924172843*phiUx[0]+623288.8754085057*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(151200.0*rdx2SqVolSq[0]*phiUy[1]+606900.0*rdx2SqVolSq[0]*phiUx[1]+((-1116705.117163882*phiUy[0])-615692.1005665087*phiUx[0]+1761980.645523667*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(20295.0*rdx2SqVolCu[0]*phiUy[1]+168555.0*rdx2SqVolCu[0]*phiUx[1]+((-110142.8429041125*phiUy[0])-170220.7572154465*phiUx[0]+522469.6620015367*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+((-204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1])-143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0])*rho[2]+(623288.8754085057*rdx2SqVolCu[1]+270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]+31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiUx[3]+(278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1]+182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]+25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1]+812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]+176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+((-244872.0*rdx2SqVol[0]*rdx2SqVolCu[1])+383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1099921.544838539*phiUy[1]*rdx2SqVolR4[1]+(1661713.956324312*rdx2SqVol[0]*phiUy[1]-376570.3622259722*rdx2SqVol[0]*phiUx[1]+((-304893.0*phiUy[0])+427329.0*phiUx[0]+839664.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(453260.3758326994*rdx2SqVolSq[0]*phiUy[1]+15016.88050162216*rdx2SqVolSq[0]*phiUx[1]+(693600.0*bcVals[0]-191520.0*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(21091.18268376621*rdx2SqVolCu[0]*phiUy[1]+39956.68007980643*rdx2SqVolCu[0]*phiUx[1]+((-25707.0*phiUy[0])-39729.0*phiUx[0]+150416.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+40344.0*phiUy[0]*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiUx[1]+(293355.0*phiUy[0]+52029.0*phiUx[0]-92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+(52029.0*phiUy[0]+293355.0*phiUx[0]-828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+(40344.0*phiUx[0]-115456.0*bcVals[0])*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-121032.0*phiUy[1]*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiUx[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiUx[0]-373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiUx[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiUx[0]-1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiUx[1]-136347.039571822*bcVals[0]*rdx2SqVolCu[0]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiUx[1]+((-36430.22463559618*phiUy[0])+36430.22463559618*phiUx[0]-64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiUx[1]+((-190246.7286525579*phiUy[0])+190246.7286525579*phiUx[0]-583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUy[2]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiUx[1]+(44403.0*phiUy[0]-44403.0*phiUx[0]-87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiUx[1]+((-44403.0*phiUy[0])+44403.0*phiUx[0]-168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-624456.0*phiLy[0]*rdx2SqVolCu[1]+(100062.3072040616*rdx2SqVol[0]*phiUx[1]-1544603.07302135*rdx2SqVol[0]*phiLy[1]+((-173313.0*phiUx[0])-3368931.0*phiLy[0]-4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2951378.203030407*rdx2SqVolSq[0]*phiUx[1]-71036.59977082233*rdx2SqVolSq[0]*phiLy[1]+((-3368931.0*phiUx[0])-173313.0*phiLy[0]-1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+((-624456.0*phiUx[0])-1862784.0*bcVals[0])*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-624456.0*phiLy[1]*rdx2SqVolCu[1]+(1097649.0*rdx2SqVol[0]*phiUx[1]-722367.0*rdx2SqVol[0]*phiLy[1]+((-1100685.379244677*phiUx[0])-847040.3948826761*phiLy[0]+5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2182239.0*rdx2SqVolSq[0]*phiUx[1]-51597.0*rdx2SqVolSq[0]*phiLy[1]+((-2275410.734360502*phiUx[0])-38955.5547130316*phiLy[0]+5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiUx[1]+(733281.0298923596*bcVals[0]-366640.5149461798*phiUx[0])*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+(22491.0*rdx2SqVol[0]*phiUx[1]-497457.0*rdx2SqVol[0]*phiLy[1]+((-38955.5547130316*phiUx[0])-2275410.734360502*phiLy[0]-934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(757809.0*rdx2SqVolSq[0]*phiUx[1]-451143.0*rdx2SqVolSq[0]*phiLy[1]+((-847040.3948826761*phiUx[0])-1100685.379244677*phiLy[0]-1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[3]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+(82239.50439417786*rdx2SqVol[0]*phiUx[1]-474351.5585164657*rdx2SqVol[0]*phiLy[1]+((-82467.0*phiUx[0])-90933.0*phiLy[0]+419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(82946.18112366594*rdx2SqVolSq[0]*phiUx[1]-109228.3200777161*rdx2SqVolSq[0]*phiLy[1]+((-90933.0*phiUx[0])-82467.0*phiLy[0]-73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(17043.37994647775*rdx2SqVolCu[1]+381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]+973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1]+2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[1]+92496.0*rho[0]*rdx2SqVolCu[1]+2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0]*rho[0])*volFac+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(115456.0*rdx2SqVolR4[1]+2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]+7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1]+257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]+747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(39767.88654178142*rdx2SqVolR4[1]+930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]+2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+40344.0*phiLy[0]*rdx2SqVolR4[1]+((-10012.98571855568*rdx2SqVol[0]*phiUx[1])+310402.557275226*rdx2SqVol[0]*phiLy[1]+(17343.0*phiUx[0]+945501.0*phiLy[0]+416232.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1113691.348758712*rdx2SqVolSq[0]*phiUx[1])+1122732.653974222*rdx2SqVolSq[0]*phiLy[1]+(1286154.0*phiUx[0]+2972058.0*phiLy[0]+5155056.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-3757176.73613406*rdx2SqVolCu[0]*phiUx[1])+639329.3979374008*rdx2SqVolCu[0]*phiLy[1]+(4278411.0*phiUx[0]+1559817.0*phiLy[0]+1.3513464e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0]*phiUx[1]+(1873368.0*phiUx[0]+5588352.0*bcVals[0])*rdx2SqVolR4[0])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = (((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(85680.0*rdx2SqVol[0]*rdx2SqVolSq[1]+260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(92496.0*rdx2SqVolCu[1]+873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[1]+407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0]*rho[0])*volFac+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1]+1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]+623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(20295.0*rdx2SqVol[0]*rdx2SqVolCu[1]+151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(168555.0*rdx2SqVol[0]*rdx2SqVolCu[1]+606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+40344.0*phiLy[1]*rdx2SqVolR4[1]+((-109839.0*rdx2SqVol[0]*phiUx[1])+413649.0*rdx2SqVol[0]*phiLy[1]+(110142.8429041125*phiUx[0]+170220.7572154465*phiLy[0]-560727.200239118*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1081578.0*rdx2SqVolSq[0]*phiUx[1])+1048338.0*rdx2SqVolSq[0]*phiLy[1]+(1116705.117163882*phiUx[0]+615692.1005665087*phiLy[0]-3464794.435460782*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2599263.0*rdx2SqVolCu[0]*phiUx[1])+464373.0*rdx2SqVolCu[0]*phiLy[1]+(2717894.290068507*phiUx[0]+350599.9924172843*phiLy[0]-6136988.564971584*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0]*phiUx[1]+(1099921.544838539*phiUx[0]-2199843.089677078*bcVals[0])*rdx2SqVolR4[0])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = (((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+(17712.0*rdx2SqVolCu[1]+472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0])*rho[2]+(197904.0*rdx2SqVol[0]*rdx2SqVolSq[1]+600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+21588.28126553848*rho[0]*rdx2SqVolCu[1]+482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiUx[3]+(72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(45449.01319060734*rdx2SqVolR4[1]+1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]+4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]+3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(3321.0*rdx2SqVol[0]*rdx2SqVolCu[1]+156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]+1873368.0*rdx2SqVolR4[0])*phiUx[2]+(984.0*rdx2SqVolR4[1]-24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]-659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-2337.0*rdx2SqVol[0]*phiUx[1])+72447.0*rdx2SqVol[0]*phiLy[1]+(4047.802737288466*phiUx[0]-52621.43558475006*phiLy[0]+97147.26569492316*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((326193.6644878315*phiUx[0]-812719.8081306986*phiLy[0]+973052.2872857346*bcVals[0])*rdx2SqVolSq[0]-287280.0*rdx2SqVolSq[0]*phiUx[1])*rdx2SqVolSq[1]+((-846963.0*rdx2SqVolCu[0]*phiUx[1])-779247.0*rdx2SqVolCu[0]*phiLy[1]+(946692.2060453439*phiUx[0]-1901183.83687717*phiLy[0]+1908983.261663653*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+(31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1]+270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0])*rho[2]+(21588.28126553848*rdx2SqVolCu[1]+143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]+204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiUx[3]+(984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(150416.0*rdx2SqVol[0]*rdx2SqVolCu[1]+693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1]+453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]+1099921.544838539*rdx2SqVolR4[0])*phiUx[2]+(39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1]+15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]-376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-25636.08400282695*rdx2SqVol[0]*phiUx[1])-176754.0528615963*rdx2SqVol[0]*phiLy[1]+(25707.0*phiUx[0]+39729.0*phiLy[0]-130872.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-182447.3038660752*rdx2SqVolSq[0]*phiUx[1])-812719.8081306986*rdx2SqVolSq[0]*phiLy[1]+(191520.0*phiUx[0]-383040.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-278113.666120527*rdx2SqVolCu[0]*phiUx[1])-566001.2949481651*rdx2SqVolCu[0]*phiLy[1]+(304893.0*phiUx[0]-427329.0*phiLy[0]+244872.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-1136585.596333157*rdx2SqVolCu[1])-2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]-742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1])-381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]-17043.37994647775*rdx2SqVolCu[0])*rho[1]+2773008.0*rho[0]*rdx2SqVolCu[1]+6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0]*rho[0])*volFac+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(5588352.0*rdx2SqVolR4[1]+1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1])-1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]-310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(1649882.317257809*rdx2SqVolR4[1]+3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]+1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]+10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+1873368.0*phiLy[0]*rdx2SqVolR4[1]+((-1500934.608060923*rdx2SqVol[0]*phiUx[1])-747388.5837200081*rdx2SqVol[0]*phiLy[1]+(1559817.0*phiUx[0]+4278411.0*phiLy[0]-2773008.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-2913967.637637727*rdx2SqVolSq[0]*phiUx[1])-257521.3140693406*rdx2SqVolSq[0]*phiLy[1]+(2972058.0*phiUx[0]+1286154.0*phiLy[0]-7782592.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-930985.9693223094*rdx2SqVolCu[0]*phiUx[1])-3195.633739964578*rdx2SqVolCu[0]*phiLy[1]+(945501.0*phiUx[0]+17343.0*phiLy[0]-2659024.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolR4[0]*phiUx[1]+(40344.0*phiUx[0]-115456.0*bcVals[0])*rdx2SqVolR4[0])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-600780.0*rdx2SqVol[0]*rdx2SqVolSq[1])-197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-2773008.0*rdx2SqVolCu[1])-2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]-472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]-17712.0*rdx2SqVolCu[0])*rho[1]+1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0]*rho[0])*volFac+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1]+973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]+97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(779247.0*rdx2SqVol[0]*rdx2SqVolCu[1]-72447.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(846963.0*rdx2SqVol[0]*rdx2SqVolCu[1]+287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-1873368.0*phiLy[1]*rdx2SqVolR4[1]+(1675359.0*rdx2SqVol[0]*phiUx[1]-998973.0*rdx2SqVol[0]*phiLy[1]+((-1901183.83687717*phiUx[0])+946692.2060453439*phiLy[0]-3735659.468951633*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(659958.0*rdx2SqVolSq[0]*phiUx[1]-156186.0*rdx2SqVolSq[0]*phiLy[1]+((-812719.8081306986*phiUx[0])+326193.6644878315*phiLy[0]-4193890.830602055*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(24363.0*rdx2SqVolCu[0]*phiUx[1]-3321.0*rdx2SqVolCu[0]*phiLy[1]+((-52621.43558475006*phiUx[0])+4047.802737288466*phiLy[0]-1119902.482954654*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-984.0*rdx2SqVolR4[0]*phiUx[1]-45449.01319060734*bcVals[0]*rdx2SqVolR4[0]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+((-825552.0*rdx2SqVolCu[1])-2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]-873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0])*rho[2]+((-260100.0*rdx2SqVol[0]*rdx2SqVolSq[1])-85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+623288.8754085057*rho[0]*rdx2SqVolCu[1]+1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiUx[3]+((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-2199843.089677078*rdx2SqVolR4[1])-6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]-3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]-560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-464373.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiUx[2]+(1047816.0*rdx2SqVolR4[1]+2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+1099921.544838539*phiLy[0]*rdx2SqVolR4[1]+((-337365.0*rdx2SqVol[0]*phiUx[1])-240705.0*rdx2SqVol[0]*phiLy[1]+(350599.9924172843*phiUx[0]+2717894.290068507*phiLy[0]-623288.8754085057*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-606900.0*rdx2SqVolSq[0]*phiUx[1])-151200.0*rdx2SqVolSq[0]*phiLy[1]+(615692.1005665087*phiUx[0]+1116705.117163882*phiLy[0]-1761980.645523667*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-168555.0*rdx2SqVolCu[0]*phiUx[1])-20295.0*rdx2SqVolCu[0]*phiLy[1]+(170220.7572154465*phiUx[0]+110142.8429041125*phiLy[0]-522469.6620015367*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+((-204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1])-143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0])*rho[2]+((-623288.8754085057*rdx2SqVolCu[1])-270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]-31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiUx[3]+((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(244872.0*rdx2SqVol[0]*rdx2SqVolCu[1]-383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1]+812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]+176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1]+182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]+25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-1099921.544838539*phiLy[1]*rdx2SqVolR4[1]+(376570.3622259722*rdx2SqVol[0]*phiUx[1]-1661713.956324312*rdx2SqVol[0]*phiLy[1]+((-427329.0*phiUx[0])+304893.0*phiLy[0]-839664.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-15016.88050162216*rdx2SqVolSq[0]*phiUx[1])-453260.3758326994*rdx2SqVolSq[0]*phiLy[1]+(191520.0*phiLy[0]-693600.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-39956.68007980643*rdx2SqVolCu[0]*phiUx[1])-21091.18268376621*rdx2SqVolCu[0]*phiLy[1]+(39729.0*phiUx[0]+25707.0*phiLy[0]-150416.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-40344.0*phiLy[0]*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiUx[1]+50064.92859277839*rdx2SqVol[0]*phiLy[1]+((-52029.0*phiUx[0])-293355.0*phiLy[0]+92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+((-293355.0*phiUx[0])-52029.0*phiLy[0]+828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+(115456.0*bcVals[0]-40344.0*phiUx[0])*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+121032.0*phiLy[1]*rdx2SqVolCu[1]+((-167649.0*rdx2SqVol[0]*phiUx[1])+221031.0*rdx2SqVol[0]*phiLy[1]+(190246.7286525579*phiUx[0]-190246.7286525579*phiLy[0]+373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-11367.0*rdx2SqVolSq[0]*phiUx[1])+29889.0*rdx2SqVolSq[0]*phiLy[1]+(36430.22463559618*phiUx[0]-36430.22463559618*phiLy[0]+1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiUx[1]+136347.039571822*bcVals[0]*rdx2SqVolCu[0])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVol[0]*phiUx[1]+35055.0*rdx2SqVol[0]*phiLy[1]+((-36430.22463559618*phiUx[0])+36430.22463559618*phiLy[0]+64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(188385.0*rdx2SqVolSq[0]*phiUx[1]-35055.0*rdx2SqVolSq[0]*phiLy[1]+((-190246.7286525579*phiUx[0])+190246.7286525579*phiLy[0]+583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLy[2]+((-39128.75979378851*rdx2SqVol[0]*phiUx[1])-190246.7286525579*rdx2SqVol[0]*phiLy[1]+(44403.0*phiUx[0]-44403.0*phiLy[0]+87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(44657.46597154836*rdx2SqVolSq[0]*phiUx[1]-36430.22463559618*rdx2SqVolSq[0]*phiLy[1]+((-44403.0*phiUx[0])+44403.0*phiLy[0]+168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-624456.0*phiUy[0]*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiLx[1]-4159512.0*rdx2SqVol[0]*bcVals[1]+((-3368931.0*phiUy[0])-173313.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiLx[1]-1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+((-173313.0*phiUy[0])-3368931.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiLx[1]-1862784.0*rdx2SqVolCu[0]*bcVals[1]-624456.0*phiLx[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+624456.0*phiUy[1]*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiLx[1]+5603489.203427449*rdx2SqVol[0]*bcVals[1]+((-847040.3948826761*phiUy[0])-1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiLx[1]+5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+((-38955.5547130316*phiUy[0])-2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiLx[1]+733281.0298923596*rdx2SqVolCu[0]*bcVals[1]-366640.5149461798*phiLx[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiLx[1]-934933.3131127583*rdx2SqVol[0]*bcVals[1]+((-2275410.734360502*phiUy[0])-38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiLx[1]-1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+((-1100685.379244677*phiUy[0])-847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiLx[1]+419832.0*rdx2SqVol[0]*bcVals[1]+((-90933.0*phiUy[0])-82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiLx[1]-73032.0*rdx2SqVolSq[0]*bcVals[1]+((-82467.0*phiUy[0])-90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-17043.37994647775*rdx2SqVolCu[1])-381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]-973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1])-2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]-1136585.596333157*rdx2SqVolCu[0])*rho[1]+92496.0*rho[0]*rdx2SqVolCu[1]+2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0]*rho[0])*volFac+(307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-39767.88654178142*rdx2SqVolR4[1])-930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]-2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1])-257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]-747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-115456.0*rdx2SqVolR4[1])-2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]-7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+40344.0*phiUy[0]*rdx2SqVolR4[1]+((-310402.557275226*rdx2SqVol[0]*phiUy[1])+10012.98571855568*rdx2SqVol[0]*phiLx[1]+416232.0*rdx2SqVol[0]*bcVals[1]+(945501.0*phiUy[0]+17343.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1122732.653974222*rdx2SqVolSq[0]*phiUy[1])+1113691.348758712*rdx2SqVolSq[0]*phiLx[1]+5155056.0*rdx2SqVolSq[0]*bcVals[1]+(2972058.0*phiUy[0]+1286154.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-639329.3979374008*rdx2SqVolCu[0]*phiUy[1])+3757176.73613406*rdx2SqVolCu[0]*phiLx[1]+1.3513464e+7*rdx2SqVolCu[0]*bcVals[1]+(1559817.0*phiUy[0]+4278411.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1649882.317257809*rdx2SqVolR4[0]*phiLx[1]+5588352.0*rdx2SqVolR4[0]*bcVals[1]+1873368.0*phiLx[0]*rdx2SqVolR4[0])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-85680.0*rdx2SqVol[0]*rdx2SqVolSq[1])-260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-92496.0*rdx2SqVolCu[1])-873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]-825552.0*rdx2SqVolCu[0])*rho[1]+407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0]*rho[0])*volFac+(39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-168555.0*rdx2SqVol[0]*rdx2SqVolCu[1])-606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-20295.0*rdx2SqVol[0]*rdx2SqVolCu[1])-151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1])-1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]-623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-40344.0*phiUy[1]*rdx2SqVolR4[1]+((-413649.0*rdx2SqVol[0]*phiUy[1])+109839.0*rdx2SqVol[0]*phiLx[1]-560727.200239118*rdx2SqVol[0]*bcVals[1]+(170220.7572154465*phiUy[0]+110142.8429041125*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1048338.0*rdx2SqVolSq[0]*phiUy[1])+1081578.0*rdx2SqVolSq[0]*phiLx[1]-3464794.435460782*rdx2SqVolSq[0]*bcVals[1]+(615692.1005665087*phiUy[0]+1116705.117163882*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-464373.0*rdx2SqVolCu[0]*phiUy[1])+2599263.0*rdx2SqVolCu[0]*phiLx[1]-6136988.564971584*rdx2SqVolCu[0]*bcVals[1]+(350599.9924172843*phiUy[0]+2717894.290068507*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1047816.0*rdx2SqVolR4[0]*phiLx[1]-2199843.089677078*rdx2SqVolR4[0]*bcVals[1]+1099921.544838539*phiLx[0]*rdx2SqVolR4[0]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+((-17712.0*rdx2SqVolCu[1])-472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0])*rho[2]+((-197904.0*rdx2SqVol[0]*rdx2SqVolSq[1])-600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+21588.28126553848*rho[0]*rdx2SqVolCu[1]+482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiLx[3]+((-984.0*rdx2SqVolR4[1])+24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]+659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-3321.0*rdx2SqVol[0]*rdx2SqVolCu[1])-156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1873368.0*rdx2SqVolR4[0])*phiLx[2]+((-45449.01319060734*rdx2SqVolR4[1])-1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]-4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-72447.0*rdx2SqVol[0]*phiUy[1])+2337.0*rdx2SqVol[0]*phiLx[1]+97147.26569492316*rdx2SqVol[0]*bcVals[1]+(4047.802737288466*phiLx[0]-52621.43558475006*phiUy[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(287280.0*rdx2SqVolSq[0]*phiLx[1]+973052.2872857346*rdx2SqVolSq[0]*bcVals[1]+(326193.6644878315*phiLx[0]-812719.8081306986*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(779247.0*rdx2SqVolCu[0]*phiUy[1]+846963.0*rdx2SqVolCu[0]*phiLx[1]+1908983.261663653*rdx2SqVolCu[0]*bcVals[1]+(946692.2060453439*phiLx[0]-1901183.83687717*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+((-31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1])-270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0])*rho[2]+((-21588.28126553848*rdx2SqVolCu[1])-143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]-204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiLx[3]+((-39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1])-15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]+376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1])-453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]-1099921.544838539*rdx2SqVolR4[0])*phiLx[2]+((-150416.0*rdx2SqVol[0]*rdx2SqVolCu[1])-693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(176754.0528615963*rdx2SqVol[0]*phiUy[1]+25636.08400282695*rdx2SqVol[0]*phiLx[1]-130872.0*rdx2SqVol[0]*bcVals[1]+(39729.0*phiUy[0]+25707.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(812719.8081306986*rdx2SqVolSq[0]*phiUy[1]+182447.3038660752*rdx2SqVolSq[0]*phiLx[1]-383040.0*rdx2SqVolSq[0]*bcVals[1]+191520.0*phiLx[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(566001.2949481651*rdx2SqVolCu[0]*phiUy[1]+278113.666120527*rdx2SqVolCu[0]*phiLx[1]+244872.0*rdx2SqVolCu[0]*bcVals[1]+(304893.0*phiLx[0]-427329.0*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1136585.596333157*rdx2SqVolCu[1]+2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]+742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1]+381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[1]+2773008.0*rho[0]*rdx2SqVolCu[1]+6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0]*rho[0])*volFac+((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-1649882.317257809*rdx2SqVolR4[1])-3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]-1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]-10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1]+1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]+310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(5588352.0*rdx2SqVolR4[1]+1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1873368.0*phiUy[0]*rdx2SqVolR4[1]+(747388.5837200081*rdx2SqVol[0]*phiUy[1]+1500934.608060923*rdx2SqVol[0]*phiLx[1]+2773008.0*rdx2SqVol[0]*bcVals[1]+(4278411.0*phiUy[0]+1559817.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(257521.3140693406*rdx2SqVolSq[0]*phiUy[1]+2913967.637637727*rdx2SqVolSq[0]*phiLx[1]+7782592.0*rdx2SqVolSq[0]*bcVals[1]+(1286154.0*phiUy[0]+2972058.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3195.633739964578*rdx2SqVolCu[0]*phiUy[1]+930985.9693223094*rdx2SqVolCu[0]*phiLx[1]+2659024.0*rdx2SqVolCu[0]*bcVals[1]+(17343.0*phiUy[0]+945501.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0]*phiLx[1]+115456.0*rdx2SqVolR4[0]*bcVals[1]+40344.0*phiLx[0]*rdx2SqVolR4[0])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = (((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(600780.0*rdx2SqVol[0]*rdx2SqVolSq[1]+197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(2773008.0*rdx2SqVolCu[1]+2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]+472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[1]+1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0]*rho[0])*volFac+((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-846963.0*rdx2SqVol[0]*rdx2SqVolCu[1])-287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(72447.0*rdx2SqVolCu[0]*rdx2SqVol[1]-779247.0*rdx2SqVol[0]*rdx2SqVolCu[1])*phiLx[2]+(1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1]+973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]+97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1873368.0*phiUy[1]*rdx2SqVolR4[1]+(998973.0*rdx2SqVol[0]*phiUy[1]-1675359.0*rdx2SqVol[0]*phiLx[1]+3735659.468951633*rdx2SqVol[0]*bcVals[1]+(946692.2060453439*phiUy[0]-1901183.83687717*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(156186.0*rdx2SqVolSq[0]*phiUy[1]-659958.0*rdx2SqVolSq[0]*phiLx[1]+4193890.830602055*rdx2SqVolSq[0]*bcVals[1]+(326193.6644878315*phiUy[0]-812719.8081306986*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3321.0*rdx2SqVolCu[0]*phiUy[1]-24363.0*rdx2SqVolCu[0]*phiLx[1]+1119902.482954654*rdx2SqVolCu[0]*bcVals[1]+(4047.802737288466*phiUy[0]-52621.43558475006*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+984.0*rdx2SqVolR4[0]*phiLx[1]+45449.01319060734*rdx2SqVolR4[0]*bcVals[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = (((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+(825552.0*rdx2SqVolCu[1]+2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]+873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0])*rho[2]+(260100.0*rdx2SqVol[0]*rdx2SqVolSq[1]+85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+623288.8754085057*rho[0]*rdx2SqVolCu[1]+1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiLx[3]+((-1047816.0*rdx2SqVolR4[1])-2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]-1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(464373.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiLx[2]+((-2199843.089677078*rdx2SqVolR4[1])-6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]-3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]-560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1099921.544838539*phiUy[0]*rdx2SqVolR4[1]+(240705.0*rdx2SqVol[0]*phiUy[1]+337365.0*rdx2SqVol[0]*phiLx[1]+623288.8754085057*rdx2SqVol[0]*bcVals[1]+(2717894.290068507*phiUy[0]+350599.9924172843*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(151200.0*rdx2SqVolSq[0]*phiUy[1]+606900.0*rdx2SqVolSq[0]*phiLx[1]+1761980.645523667*rdx2SqVolSq[0]*bcVals[1]+(1116705.117163882*phiUy[0]+615692.1005665087*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(20295.0*rdx2SqVolCu[0]*phiUy[1]+168555.0*rdx2SqVolCu[0]*phiLx[1]+522469.6620015367*rdx2SqVolCu[0]*bcVals[1]+(110142.8429041125*phiUy[0]+170220.7572154465*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+(204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1]+143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0])*rho[2]+(623288.8754085057*rdx2SqVolCu[1]+270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]+31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiLx[3]+((-278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1])-182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]-25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1])-812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]-176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(244872.0*rdx2SqVol[0]*rdx2SqVolCu[1]-383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1099921.544838539*phiUy[1]*rdx2SqVolR4[1]+(1661713.956324312*rdx2SqVol[0]*phiUy[1]-376570.3622259722*rdx2SqVol[0]*phiLx[1]+839664.0*rdx2SqVol[0]*bcVals[1]+(304893.0*phiUy[0]-427329.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(453260.3758326994*rdx2SqVolSq[0]*phiUy[1]+15016.88050162216*rdx2SqVolSq[0]*phiLx[1]+693600.0*rdx2SqVolSq[0]*bcVals[1]+191520.0*phiUy[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(21091.18268376621*rdx2SqVolCu[0]*phiUy[1]+39956.68007980643*rdx2SqVolCu[0]*phiLx[1]+150416.0*rdx2SqVolCu[0]*bcVals[1]+(25707.0*phiUy[0]+39729.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-40344.0*phiUy[0]*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiLx[1]-92496.0*rdx2SqVol[0]*bcVals[1]+((-293355.0*phiUy[0])-52029.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiLx[1]-828720.0*rdx2SqVolSq[0]*bcVals[1]+((-52029.0*phiUy[0])-293355.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiLx[1]-115456.0*rdx2SqVolCu[0]*bcVals[1]-40344.0*phiLx[0]*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-121032.0*phiUy[1]*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiLx[1]-373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiLx[1]-1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiLx[1]-136347.039571822*rdx2SqVolCu[0]*bcVals[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiLx[1]-64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiLx[1]-583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiLx[1]-87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiLx[1]-168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiUy[0]-44403.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+624456.0*phiLy[0]*rdx2SqVolCu[1]+((-1544603.07302135*rdx2SqVol[0]*phiLy[1])+100062.3072040616*rdx2SqVol[0]*phiLx[1]+4159512.0*rdx2SqVol[0]*bcVals[1]+(3368931.0*phiLy[0]+173313.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-71036.59977082233*rdx2SqVolSq[0]*phiLy[1])+2951378.203030407*rdx2SqVolSq[0]*phiLx[1]+1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+(173313.0*phiLy[0]+3368931.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiLx[1]+1862784.0*rdx2SqVolCu[0]*bcVals[1]+624456.0*phiLx[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]-624456.0*phiLy[1]*rdx2SqVolCu[1]+((-722367.0*rdx2SqVol[0]*phiLy[1])+1097649.0*rdx2SqVol[0]*phiLx[1]-5603489.203427449*rdx2SqVol[0]*bcVals[1]+(847040.3948826761*phiLy[0]+1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-51597.0*rdx2SqVolSq[0]*phiLy[1])+2182239.0*rdx2SqVolSq[0]*phiLx[1]-5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+(38955.5547130316*phiLy[0]+2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiLx[1]-733281.0298923596*rdx2SqVolCu[0]*bcVals[1]+366640.5149461798*phiLx[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+((-497457.0*rdx2SqVol[0]*phiLy[1])+22491.0*rdx2SqVol[0]*phiLx[1]+934933.3131127583*rdx2SqVol[0]*bcVals[1]+(2275410.734360502*phiLy[0]+38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-451143.0*rdx2SqVolSq[0]*phiLy[1])+757809.0*rdx2SqVolSq[0]*phiLx[1]+1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+(1100685.379244677*phiLy[0]+847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+((-474351.5585164657*rdx2SqVol[0]*phiLy[1])+82239.50439417786*rdx2SqVol[0]*phiLx[1]-419832.0*rdx2SqVol[0]*bcVals[1]+(90933.0*phiLy[0]+82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-109228.3200777161*rdx2SqVolSq[0]*phiLy[1])+82946.18112366594*rdx2SqVolSq[0]*phiLx[1]+73032.0*rdx2SqVolSq[0]*bcVals[1]+(82467.0*phiLy[0]+90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-17043.37994647775*rdx2SqVolCu[1])-381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]-973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1]+2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[1]-92496.0*rho[0]*rdx2SqVolCu[1]-2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0]*rho[0])*volFac+(307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-115456.0*rdx2SqVolR4[1])-2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]-7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-39767.88654178142*rdx2SqVolR4[1])-930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]-2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1])-257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]-747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]-40344.0*phiLy[0]*rdx2SqVolR4[1]+(310402.557275226*rdx2SqVol[0]*phiLy[1]-10012.98571855568*rdx2SqVol[0]*phiLx[1]-416232.0*rdx2SqVol[0]*bcVals[1]+((-945501.0*phiLy[0])-17343.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1122732.653974222*rdx2SqVolSq[0]*phiLy[1]-1113691.348758712*rdx2SqVolSq[0]*phiLx[1]-5155056.0*rdx2SqVolSq[0]*bcVals[1]+((-2972058.0*phiLy[0])-1286154.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(639329.3979374008*rdx2SqVolCu[0]*phiLy[1]-3757176.73613406*rdx2SqVolCu[0]*phiLx[1]-1.3513464e+7*rdx2SqVolCu[0]*bcVals[1]+((-1559817.0*phiLy[0])-4278411.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0]*phiLx[1]-5588352.0*rdx2SqVolR4[0]*bcVals[1]-1873368.0*phiLx[0]*rdx2SqVolR4[0]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = (((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-85680.0*rdx2SqVol[0]*rdx2SqVolSq[1])-260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(92496.0*rdx2SqVolCu[1]+873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[1]-407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0]*rho[0])*volFac+(39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1])-1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]-623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-168555.0*rdx2SqVol[0]*rdx2SqVolCu[1])-606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-20295.0*rdx2SqVol[0]*rdx2SqVolCu[1])-151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+40344.0*phiLy[1]*rdx2SqVolR4[1]+(413649.0*rdx2SqVol[0]*phiLy[1]-109839.0*rdx2SqVol[0]*phiLx[1]+560727.200239118*rdx2SqVol[0]*bcVals[1]+((-170220.7572154465*phiLy[0])-110142.8429041125*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1048338.0*rdx2SqVolSq[0]*phiLy[1]-1081578.0*rdx2SqVolSq[0]*phiLx[1]+3464794.435460782*rdx2SqVolSq[0]*bcVals[1]+((-615692.1005665087*phiLy[0])-1116705.117163882*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(464373.0*rdx2SqVolCu[0]*phiLy[1]-2599263.0*rdx2SqVolCu[0]*phiLx[1]+6136988.564971584*rdx2SqVolCu[0]*bcVals[1]+((-350599.9924172843*phiLy[0])-2717894.290068507*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0]*phiLx[1]+2199843.089677078*rdx2SqVolR4[0]*bcVals[1]-1099921.544838539*phiLx[0]*rdx2SqVolR4[0])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+((-17712.0*rdx2SqVolCu[1])-472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0])*rho[2]+(197904.0*rdx2SqVol[0]*rdx2SqVolSq[1]+600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-21588.28126553848*rho[0]*rdx2SqVolCu[1]-482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiLx[3]+((-45449.01319060734*rdx2SqVolR4[1])-1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]-4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-984.0*rdx2SqVolR4[1])+24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]+659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-3321.0*rdx2SqVol[0]*rdx2SqVolCu[1])-156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1873368.0*rdx2SqVolR4[0])*phiLx[2]+(72447.0*rdx2SqVol[0]*phiLy[1]-2337.0*rdx2SqVol[0]*phiLx[1]-97147.26569492316*rdx2SqVol[0]*bcVals[1]+(52621.43558475006*phiLy[0]-4047.802737288466*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-287280.0*rdx2SqVolSq[0]*phiLx[1])-973052.2872857346*rdx2SqVolSq[0]*bcVals[1]+(812719.8081306986*phiLy[0]-326193.6644878315*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-779247.0*rdx2SqVolCu[0]*phiLy[1])-846963.0*rdx2SqVolCu[0]*phiLx[1]-1908983.261663653*rdx2SqVolCu[0]*bcVals[1]+(1901183.83687717*phiLy[0]-946692.2060453439*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+((-31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1])-270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0])*rho[2]+(21588.28126553848*rdx2SqVolCu[1]+143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]+204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiLx[3]+((-150416.0*rdx2SqVol[0]*rdx2SqVolCu[1])-693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1])-15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]+376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1])-453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]-1099921.544838539*rdx2SqVolR4[0])*phiLx[2]+((-176754.0528615963*rdx2SqVol[0]*phiLy[1])-25636.08400282695*rdx2SqVol[0]*phiLx[1]+130872.0*rdx2SqVol[0]*bcVals[1]+((-39729.0*phiLy[0])-25707.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-812719.8081306986*rdx2SqVolSq[0]*phiLy[1])-182447.3038660752*rdx2SqVolSq[0]*phiLx[1]+383040.0*rdx2SqVolSq[0]*bcVals[1]-191520.0*phiLx[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-566001.2949481651*rdx2SqVolCu[0]*phiLy[1])-278113.666120527*rdx2SqVolCu[0]*phiLx[1]-244872.0*rdx2SqVolCu[0]*bcVals[1]+(427329.0*phiLy[0]-304893.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1136585.596333157*rdx2SqVolCu[1]+2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]+742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1])-381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]-17043.37994647775*rdx2SqVolCu[0])*rho[1]-2773008.0*rho[0]*rdx2SqVolCu[1]-6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0]*rho[0])*volFac+((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-5588352.0*rdx2SqVolR4[1])-1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1649882.317257809*rdx2SqVolR4[1])-3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]-1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]-10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1]+1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]+310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]-1873368.0*phiLy[0]*rdx2SqVolR4[1]+((-747388.5837200081*rdx2SqVol[0]*phiLy[1])-1500934.608060923*rdx2SqVol[0]*phiLx[1]-2773008.0*rdx2SqVol[0]*bcVals[1]+((-4278411.0*phiLy[0])-1559817.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-257521.3140693406*rdx2SqVolSq[0]*phiLy[1])-2913967.637637727*rdx2SqVolSq[0]*phiLx[1]-7782592.0*rdx2SqVolSq[0]*bcVals[1]+((-1286154.0*phiLy[0])-2972058.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-3195.633739964578*rdx2SqVolCu[0]*phiLy[1])-930985.9693223094*rdx2SqVolCu[0]*phiLx[1]-2659024.0*rdx2SqVolCu[0]*bcVals[1]+((-17343.0*phiLy[0])-945501.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolR4[0]*phiLx[1]-115456.0*rdx2SqVolR4[0]*bcVals[1]-40344.0*phiLx[0]*rdx2SqVolR4[0]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(600780.0*rdx2SqVol[0]*rdx2SqVolSq[1]+197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-2773008.0*rdx2SqVolCu[1])-2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]-472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]-17712.0*rdx2SqVolCu[0])*rho[1]-1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0]*rho[0])*volFac+((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1])-973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]-97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-846963.0*rdx2SqVol[0]*rdx2SqVolCu[1])-287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(72447.0*rdx2SqVolCu[0]*rdx2SqVol[1]-779247.0*rdx2SqVol[0]*rdx2SqVolCu[1])*phiLx[2]-1873368.0*phiLy[1]*rdx2SqVolR4[1]+((-998973.0*rdx2SqVol[0]*phiLy[1])+1675359.0*rdx2SqVol[0]*phiLx[1]-3735659.468951633*rdx2SqVol[0]*bcVals[1]+(1901183.83687717*phiLx[0]-946692.2060453439*phiLy[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-156186.0*rdx2SqVolSq[0]*phiLy[1])+659958.0*rdx2SqVolSq[0]*phiLx[1]-4193890.830602055*rdx2SqVolSq[0]*bcVals[1]+(812719.8081306986*phiLx[0]-326193.6644878315*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-3321.0*rdx2SqVolCu[0]*phiLy[1])+24363.0*rdx2SqVolCu[0]*phiLx[1]-1119902.482954654*rdx2SqVolCu[0]*bcVals[1]+(52621.43558475006*phiLx[0]-4047.802737288466*phiLy[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-984.0*rdx2SqVolR4[0]*phiLx[1]-45449.01319060734*rdx2SqVolR4[0]*bcVals[1]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = (((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+(825552.0*rdx2SqVolCu[1]+2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]+873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0])*rho[2]+((-260100.0*rdx2SqVol[0]*rdx2SqVolSq[1])-85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-623288.8754085057*rho[0]*rdx2SqVolCu[1]-1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiLx[3]+(2199843.089677078*rdx2SqVolR4[1]+6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]+3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]+560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1047816.0*rdx2SqVolR4[1])-2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]-1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(464373.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiLx[2]-1099921.544838539*phiLy[0]*rdx2SqVolR4[1]+((-240705.0*rdx2SqVol[0]*phiLy[1])-337365.0*rdx2SqVol[0]*phiLx[1]-623288.8754085057*rdx2SqVol[0]*bcVals[1]+((-2717894.290068507*phiLy[0])-350599.9924172843*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-151200.0*rdx2SqVolSq[0]*phiLy[1])-606900.0*rdx2SqVolSq[0]*phiLx[1]-1761980.645523667*rdx2SqVolSq[0]*bcVals[1]+((-1116705.117163882*phiLy[0])-615692.1005665087*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-20295.0*rdx2SqVolCu[0]*phiLy[1])-168555.0*rdx2SqVolCu[0]*phiLx[1]-522469.6620015367*rdx2SqVolCu[0]*bcVals[1]+((-110142.8429041125*phiLy[0])-170220.7572154465*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+(204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1]+143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0])*rho[2]+((-623288.8754085057*rdx2SqVolCu[1])-270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]-31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiLx[3]+((-244872.0*rdx2SqVol[0]*rdx2SqVolCu[1])+383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1])-182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]-25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1])-812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]-176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]-1099921.544838539*phiLy[1]*rdx2SqVolR4[1]+((-1661713.956324312*rdx2SqVol[0]*phiLy[1])+376570.3622259722*rdx2SqVol[0]*phiLx[1]-839664.0*rdx2SqVol[0]*bcVals[1]+(427329.0*phiLx[0]-304893.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-453260.3758326994*rdx2SqVolSq[0]*phiLy[1])-15016.88050162216*rdx2SqVolSq[0]*phiLx[1]-693600.0*rdx2SqVolSq[0]*bcVals[1]-191520.0*phiLy[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-21091.18268376621*rdx2SqVolCu[0]*phiLy[1])-39956.68007980643*rdx2SqVolCu[0]*phiLx[1]-150416.0*rdx2SqVolCu[0]*bcVals[1]+((-25707.0*phiLy[0])-39729.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+40344.0*phiLy[0]*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiLy[1]+50064.92859277839*rdx2SqVol[0]*phiLx[1]+92496.0*rdx2SqVol[0]*bcVals[1]+(293355.0*phiLy[0]+52029.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+288932.0554646022*rdx2SqVolSq[0]*phiLx[1]+828720.0*rdx2SqVolSq[0]*bcVals[1]+(52029.0*phiLy[0]+293355.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiLx[1]+115456.0*rdx2SqVolCu[0]*bcVals[1]+40344.0*phiLx[0]*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+121032.0*phiLy[1]*rdx2SqVolCu[1]+(221031.0*rdx2SqVol[0]*phiLy[1]-167649.0*rdx2SqVol[0]*phiLx[1]+373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLy[0]-190246.7286525579*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(29889.0*rdx2SqVolSq[0]*phiLy[1]-11367.0*rdx2SqVolSq[0]*phiLx[1]+1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLy[0]-36430.22463559618*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiLx[1]+136347.039571822*rdx2SqVolCu[0]*bcVals[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+(35055.0*rdx2SqVol[0]*phiLy[1]+35055.0*rdx2SqVol[0]*phiLx[1]+64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-35055.0*rdx2SqVolSq[0]*phiLy[1])+188385.0*rdx2SqVolSq[0]*phiLx[1]+583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-190246.7286525579*rdx2SqVol[0]*phiLy[1])-39128.75979378851*rdx2SqVol[0]*phiLx[1]+87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLy[0]-44403.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-36430.22463559618*rdx2SqVolSq[0]*phiLy[1])+44657.46597154836*rdx2SqVolSq[0]*phiLx[1]+168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (16.0*rho[0]*omega*volFac+((-8.660254037844386*rdx2SqVol[1]*phiUy[2])+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0])*omega+18.0*phiC[0]*rdx2SqVol[1]+18.0*phiC[0]*rdx2SqVol[0])/(18.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[1] = (16.0*rho[1]*omega*volFac+((-8.660254037844386*rdx2SqVol[1]*phiUy[3])+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1]-18.0*phiC[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])*omega+18.0*phiC[1]*rdx2SqVol[1]+46.0*rdx2SqVol[0]*phiC[1])/(18.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 
  phiC[2] = (16.0*rho[2]*omega*volFac+((-8.660254037844386*rdx2SqVol[0]*phiUx[3])+8.660254037844386*rdx2SqVol[0]*phiLx[3]-7.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]-7.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-46.0*rdx2SqVol[1])-18.0*rdx2SqVol[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1])*omega+(46.0*rdx2SqVol[1]+18.0*rdx2SqVol[0])*phiC[2])/(46.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[3] = (16.0*rho[3]*omega*volFac+((-7.0*rdx2SqVol[1]*phiUy[3])-7.0*rdx2SqVol[0]*phiUx[3]-7.0*rdx2SqVol[1]*phiLy[3]-7.0*rdx2SqVol[0]*phiLx[3]+((-46.0*rdx2SqVol[1])-46.0*rdx2SqVol[0])*phiC[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1])*omega+(46.0*rdx2SqVol[1]+46.0*rdx2SqVol[0])*phiC[3])/(46.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[0]*rho[1]+288.0*rho[0]*rdx2SqVol[1]+2096.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(162.0*phiUy[0]+162.0*phiLy[0]-324.0*phiC[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]-31.17691453623978*rdx2SqVol[0]*phiUx[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]+(1179.0*phiUy[0]+54.0*phiUx[0]+1179.0*phiLy[0]-3060.0*phiC[0]+1296.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiUx[1]+(1416.0*phiUx[0]-3528.0*phiC[0]+4224.0*bcVals[0])*rdx2SqVolSq[0])*omega+324.0*phiC[0]*rdx2SqVolSq[1]+3060.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*phiC[0]*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]+471.1178196587346*rdx2SqVol[0]*rho[0])*omega*volFac+(((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]+162.0*phiLy[1]-324.0*phiC[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]-342.0*rdx2SqVol[0]*phiUx[1]+351.0*rdx2SqVol[0]*phiLy[1]-3060.0*rdx2SqVol[0]*phiC[1]+(265.0037735580381*phiUy[0]+342.9460598986376*phiUx[0]+265.0037735580381*phiLy[0]-1745.907214029428*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiUx[1]-3528.0*rdx2SqVolSq[0]*phiC[1]+(831.384387633061*phiUx[0]-1662.768775266122*bcVals[0])*rdx2SqVolSq[0])*omega+324.0*phiC[1]*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]*phiC[1])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = ((859.0972005541631*rdx2SqVol[0]*rho[3]+(736.0*rdx2SqVol[1]+2096.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(138.0*rdx2SqVol[0]*rdx2SqVol[1]+1416.0*rdx2SqVolSq[0])*phiUx[2]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiC[2]+(398.3716857408418*phiUy[0]-398.3716857408418*phiLy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiUy[0]-1134.493278957615*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiC[2])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[0]*rho[2])*omega*volFac+(((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiUx[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiC[3]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+(876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1]+831.384387633061*rdx2SqVolSq[0])*phiUx[2]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiUy[0]-255.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiC[3])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[1]-864.0*rho[0]*rdx2SqVol[1]-2256.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-486.0*phiUy[0])-486.0*phiLy[0]+972.0*phiC[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+467.6537180435967*rdx2SqVol[0]*phiUx[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+((-1269.0*phiUy[0])-486.0*phiUx[0]-1269.0*phiLy[0]+3024.0*phiC[0]+864.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiUx[1]+((-984.0*phiUx[0])+984.0*phiC[0]+2816.0*bcVals[0])*rdx2SqVolSq[0])*omega-972.0*phiC[0]*rdx2SqVolSq[1]-3024.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-984.0*phiC[0]*rdx2SqVolSq[0]))/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]-526.5434455009387*rdx2SqVol[0]*rho[0])*omega*volFac+(((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]+486.0*phiLy[1]-972.0*phiC[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]-522.0*rdx2SqVol[0]*phiUx[1]+243.0*rdx2SqVol[0]*phiLy[1]-3024.0*rdx2SqVol[0]*phiC[1]+((-296.1806880942779*phiUy[0])+592.3613761885558*phiUx[0]-296.1806880942779*phiLy[0]+1163.938142686285*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiUx[1]-984.0*rdx2SqVolSq[0]*phiC[1]+1108.512516844081*bcVals[0]*rdx2SqVolSq[0])*omega+972.0*phiC[1]*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]*phiC[1])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[3]+((-2208.0*rdx2SqVol[1])-2256.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-1242.0*rdx2SqVol[0]*rdx2SqVol[1])-984.0*rdx2SqVolSq[0])*phiUx[2]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiC[2]+(1195.115057222525*phiLy[0]-1195.115057222525*phiUy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiLy[0]-1221.095819336058*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiC[2]))/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[0]*rho[2])*omega*volFac+(((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUx[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiC[3]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiLy[0]-285.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiC[3])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[1]-288.0*rho[0]*rdx2SqVol[1]-2096.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-162.0*phiUy[0])-162.0*phiLy[0]+324.0*phiC[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]-31.17691453623978*rdx2SqVol[0]*phiLx[1]-1296.0*rdx2SqVol[0]*bcVals[1]+((-1179.0*phiUy[0])-1179.0*phiLy[0]-54.0*phiLx[0]+3060.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiLx[1]-4224.0*rdx2SqVolSq[0]*bcVals[1]+(3528.0*phiC[0]-1416.0*phiLx[0])*rdx2SqVolSq[0])*omega-324.0*phiC[0]*rdx2SqVolSq[1]-3060.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*phiC[0]*rdx2SqVolSq[0]))/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]-471.1178196587346*rdx2SqVol[0]*rho[0])*omega*volFac+(((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]+162.0*phiLy[1]-324.0*phiC[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]+351.0*rdx2SqVol[0]*phiLy[1]-342.0*rdx2SqVol[0]*phiLx[1]-3060.0*rdx2SqVol[0]*phiC[1]+1745.907214029428*rdx2SqVol[0]*bcVals[1]+((-265.0037735580381*phiUy[0])-265.0037735580381*phiLy[0]-342.9460598986376*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiLx[1]-3528.0*rdx2SqVolSq[0]*phiC[1]+1662.768775266122*rdx2SqVolSq[0]*bcVals[1]-831.384387633061*phiLx[0]*rdx2SqVolSq[0])*omega+324.0*phiC[1]*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]*phiC[1])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[3]+((-736.0*rdx2SqVol[1])-2096.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiLx[3]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-138.0*rdx2SqVol[0]*rdx2SqVol[1])-1416.0*rdx2SqVolSq[0])*phiLx[2]+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiC[2]+(398.3716857408418*phiLy[0]-398.3716857408418*phiUy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiLy[0]-1134.493278957615*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiC[2]))/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[0]*rho[2])*omega*volFac+(((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiLx[3]+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiC[3]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+((-876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1])-831.384387633061*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiLy[0]-255.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiC[3])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[0]*rho[1]+864.0*rho[0]*rdx2SqVol[1]+2256.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(486.0*phiUy[0]+486.0*phiLy[0]-972.0*phiC[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+467.6537180435967*rdx2SqVol[0]*phiLx[1]+864.0*rdx2SqVol[0]*bcVals[1]+(1269.0*phiUy[0]+1269.0*phiLy[0]+486.0*phiLx[0]-3024.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiLx[1]+2816.0*rdx2SqVolSq[0]*bcVals[1]+(984.0*phiLx[0]-984.0*phiC[0])*rdx2SqVolSq[0])*omega+972.0*phiC[0]*rdx2SqVolSq[1]+3024.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+984.0*phiC[0]*rdx2SqVolSq[0])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]+526.5434455009387*rdx2SqVol[0]*rho[0])*omega*volFac+(((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]+486.0*phiLy[1]-972.0*phiC[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]+243.0*rdx2SqVol[0]*phiLy[1]-522.0*rdx2SqVol[0]*phiLx[1]-3024.0*rdx2SqVol[0]*phiC[1]+1163.938142686285*rdx2SqVol[0]*bcVals[1]+(296.1806880942779*phiUy[0]+296.1806880942779*phiLy[0]-592.3613761885558*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiLx[1]-984.0*rdx2SqVolSq[0]*phiC[1]+1108.512516844081*rdx2SqVolSq[0]*bcVals[1])*omega+972.0*phiC[1]*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]*phiC[1])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = ((415.6921938165305*rdx2SqVol[0]*rho[3]+(2208.0*rdx2SqVol[1]+2256.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiLx[3]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1242.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiLx[2]+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiC[2]+(1195.115057222525*phiUy[0]-1195.115057222525*phiLy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiUy[0]-1221.095819336058*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiC[2])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[0]*rho[2])*omega*volFac+(((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLx[3]+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiC[3]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]-1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiUy[0]-285.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiC[3])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[1]*rho[2]+2096.0*rho[0]*rdx2SqVol[1]+288.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(4224.0*rdx2SqVolSq[1]+1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(1416.0*phiUy[0]-3528.0*phiC[0])*rdx2SqVolSq[1]+((-1134.493278957615*rdx2SqVol[0]*phiUx[1])+1134.493278957615*rdx2SqVol[0]*phiLx[1]+(54.0*phiUy[0]+1179.0*phiUx[0]+1179.0*phiLx[0]-3060.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-155.8845726811989*rdx2SqVolSq[0]*phiUx[1]+155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+(162.0*phiUx[0]+162.0*phiLx[0]-324.0*phiC[0])*rdx2SqVolSq[0])*omega+3528.0*phiC[0]*rdx2SqVolSq[1]+3060.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+324.0*phiC[0]*rdx2SqVolSq[0])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = ((859.0972005541631*rdx2SqVol[1]*rho[3]+(2096.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[1])*omega*volFac+(((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(1416.0*phiUy[1]-3528.0*phiC[1])*rdx2SqVolSq[1]+(138.0*rdx2SqVol[0]*phiUy[1]-917.0*rdx2SqVol[0]*phiUx[1]-917.0*rdx2SqVol[0]*phiLx[1]-7820.0*rdx2SqVol[0]*phiC[1]+(1134.493278957615*phiUx[0]-1134.493278957615*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-322.0*rdx2SqVolSq[0]*phiUx[1]-322.0*rdx2SqVolSq[0]*phiLx[1]-2116.0*rdx2SqVolSq[0]*phiC[1]+(398.3716857408418*phiUx[0]-398.3716857408418*phiLx[0])*rdx2SqVolSq[0])*omega+3528.0*phiC[1]*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]*phiC[1])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]+471.1178196587346*rho[0]*rdx2SqVol[1])*omega*volFac+(((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]+((-3528.0*rdx2SqVolSq[1])-3060.0*rdx2SqVol[0]*rdx2SqVol[1]-324.0*rdx2SqVolSq[0])*phiC[2]+((-1662.768775266122*rdx2SqVolSq[1])-1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+831.384387633061*phiUy[0]*rdx2SqVolSq[1]+((-255.0*rdx2SqVol[0]*phiUx[1])+255.0*rdx2SqVol[0]*phiLx[1]+(342.9460598986376*phiUy[0]+265.0037735580381*phiUx[0]+265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0])*phiC[2])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[1]*rho[1])*omega*volFac+(((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+((-3528.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-2116.0*rdx2SqVolSq[0])*phiC[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]+831.384387633061*phiUy[1]*rdx2SqVolSq[1]+(876.4177086298519*rdx2SqVol[0]*phiUy[1]-206.1140461006964*rdx2SqVol[0]*phiUx[1]-206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiUx[0]-255.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0])*phiC[3])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[2]-2256.0*rho[0]*rdx2SqVol[1]-864.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(984.0*phiC[0]-984.0*phiUy[0])*rdx2SqVolSq[1]+(1221.095819336058*rdx2SqVol[0]*phiUx[1]-1221.095819336058*rdx2SqVol[0]*phiLx[1]+((-486.0*phiUy[0])-1269.0*phiUx[0]-1269.0*phiLx[0]+3024.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0]*phiUx[1]-467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+((-486.0*phiUx[0])-486.0*phiLx[0]+972.0*phiC[0])*rdx2SqVolSq[0])*omega-984.0*phiC[0]*rdx2SqVolSq[1]-3024.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-972.0*phiC[0]*rdx2SqVolSq[0]))/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[3]+((-2256.0*rdx2SqVol[1])-2208.0*rdx2SqVol[0])*rho[1])*omega*volFac+((969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(984.0*phiC[1]-984.0*phiUy[1])*rdx2SqVolSq[1]+((-1242.0*rdx2SqVol[0]*phiUy[1])+987.0*rdx2SqVol[0]*phiUx[1]+987.0*rdx2SqVol[0]*phiLx[1]+7728.0*rdx2SqVol[0]*phiC[1]+(1221.095819336058*phiLx[0]-1221.095819336058*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+966.0*rdx2SqVolSq[0]*phiUx[1]+966.0*rdx2SqVolSq[0]*phiLx[1]+6348.0*rdx2SqVolSq[0]*phiC[1]+(1195.115057222525*phiLx[0]-1195.115057222525*phiUx[0])*rdx2SqVolSq[0])*omega-984.0*phiC[1]*rdx2SqVolSq[1]-7728.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]-6348.0*rdx2SqVolSq[0]*phiC[1]))/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]-526.5434455009387*rho[0]*rdx2SqVol[1])*omega*volFac+(((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+((-984.0*rdx2SqVolSq[1])-3024.0*rdx2SqVol[0]*rdx2SqVol[1]-972.0*rdx2SqVolSq[0])*phiC[2]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(285.0*rdx2SqVol[0]*phiUx[1]-285.0*rdx2SqVol[0]*phiLx[1]+(592.3613761885558*phiUy[0]-296.1806880942779*phiUx[0]-296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0])*phiC[2])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[1]*rho[1])*omega*volFac+((24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+((-984.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-6348.0*rdx2SqVolSq[0])*phiC[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+(1513.812405815199*rdx2SqVol[0]*phiUy[1]+230.3627574066607*rdx2SqVol[0]*phiUx[1]+230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiLx[0]-285.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0])*phiC[3])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[2]-2096.0*rho[0]*rdx2SqVol[1]-288.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-4224.0*rdx2SqVolSq[1])-1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(3528.0*phiC[0]-1416.0*phiLy[0])*rdx2SqVolSq[1]+(1134.493278957615*rdx2SqVol[0]*phiUx[1]-1134.493278957615*rdx2SqVol[0]*phiLx[1]+((-1179.0*phiUx[0])-54.0*phiLy[0]-1179.0*phiLx[0]+3060.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0]*phiUx[1]-155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+((-162.0*phiUx[0])-162.0*phiLx[0]+324.0*phiC[0])*rdx2SqVolSq[0])*omega-3528.0*phiC[0]*rdx2SqVolSq[1]-3060.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-324.0*phiC[0]*rdx2SqVolSq[0]))/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[3]+((-2096.0*rdx2SqVol[1])-736.0*rdx2SqVol[0])*rho[1])*omega*volFac+((-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(3528.0*phiC[1]-1416.0*phiLy[1])*rdx2SqVolSq[1]+(917.0*rdx2SqVol[0]*phiUx[1]-138.0*rdx2SqVol[0]*phiLy[1]+917.0*rdx2SqVol[0]*phiLx[1]+7820.0*rdx2SqVol[0]*phiC[1]+(1134.493278957615*phiLx[0]-1134.493278957615*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+322.0*rdx2SqVolSq[0]*phiUx[1]+322.0*rdx2SqVolSq[0]*phiLx[1]+2116.0*rdx2SqVolSq[0]*phiC[1]+(398.3716857408418*phiLx[0]-398.3716857408418*phiUx[0])*rdx2SqVolSq[0])*omega-3528.0*phiC[1]*rdx2SqVolSq[1]-7820.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]-2116.0*rdx2SqVolSq[0]*phiC[1]))/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]-471.1178196587346*rho[0]*rdx2SqVol[1])*omega*volFac+(((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+(1662.768775266122*rdx2SqVolSq[1]+1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]+((-3528.0*rdx2SqVolSq[1])-3060.0*rdx2SqVol[0]*rdx2SqVol[1]-324.0*rdx2SqVolSq[0])*phiC[2]-831.384387633061*phiLy[0]*rdx2SqVolSq[1]+(255.0*rdx2SqVol[0]*phiUx[1]-255.0*rdx2SqVol[0]*phiLx[1]+((-265.0037735580381*phiUx[0])-342.9460598986376*phiLy[0]-265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0])*phiC[2])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[1]*rho[1])*omega*volFac+(((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+((-3528.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-2116.0*rdx2SqVolSq[0])*phiC[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[1]*rdx2SqVolSq[1]+(206.1140461006964*rdx2SqVol[0]*phiUx[1]-876.4177086298519*rdx2SqVol[0]*phiLy[1]+206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiLx[0]-255.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0])*phiC[3])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[1]*rho[2]+2256.0*rho[0]*rdx2SqVol[1]+864.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(984.0*phiLy[0]-984.0*phiC[0])*rdx2SqVolSq[1]+((-1221.095819336058*rdx2SqVol[0]*phiUx[1])+1221.095819336058*rdx2SqVol[0]*phiLx[1]+(1269.0*phiUx[0]+486.0*phiLy[0]+1269.0*phiLx[0]-3024.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-467.6537180435967*rdx2SqVolSq[0]*phiUx[1]+467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+(486.0*phiUx[0]+486.0*phiLx[0]-972.0*phiC[0])*rdx2SqVolSq[0])*omega+984.0*phiC[0]*rdx2SqVolSq[1]+3024.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+972.0*phiC[0]*rdx2SqVolSq[0])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = ((415.6921938165305*rdx2SqVol[1]*rho[3]+(2256.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[1])*omega*volFac+((-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+(969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(984.0*phiLy[1]-984.0*phiC[1])*rdx2SqVolSq[1]+((-987.0*rdx2SqVol[0]*phiUx[1])+1242.0*rdx2SqVol[0]*phiLy[1]-987.0*rdx2SqVol[0]*phiLx[1]-7728.0*rdx2SqVol[0]*phiC[1]+(1221.095819336058*phiUx[0]-1221.095819336058*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-966.0*rdx2SqVolSq[0]*phiUx[1]-966.0*rdx2SqVolSq[0]*phiLx[1]-6348.0*rdx2SqVolSq[0]*phiC[1]+(1195.115057222525*phiUx[0]-1195.115057222525*phiLx[0])*rdx2SqVolSq[0])*omega+984.0*phiC[1]*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]*phiC[1])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]+526.5434455009387*rho[0]*rdx2SqVol[1])*omega*volFac+(((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+((-984.0*rdx2SqVolSq[1])-3024.0*rdx2SqVol[0]*rdx2SqVol[1]-972.0*rdx2SqVolSq[0])*phiC[2]+((-285.0*rdx2SqVol[0]*phiUx[1])+285.0*rdx2SqVol[0]*phiLx[1]+(296.1806880942779*phiUx[0]-592.3613761885558*phiLy[0]+296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0])*phiC[2])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[1]*rho[1])*omega*volFac+(((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+(24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+((-984.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-6348.0*rdx2SqVolSq[0])*phiC[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+((-230.3627574066607*rdx2SqVol[0]*phiUx[1])-1513.812405815199*rdx2SqVol[0]*phiLy[1]-230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiUx[0]-285.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0])*phiC[3])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(624456.0*phiUy[0]-1555848.0*phiC[0])*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiUx[1]+(3368931.0*phiUy[0]+173313.0*phiUx[0]-1.1189052e+7*phiC[0]+4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiUx[1]+(173313.0*phiUy[0]+3368931.0*phiUx[0]-1.1189052e+7*phiC[0]+1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+(624456.0*phiUx[0]-1555848.0*phiC[0]+1862784.0*bcVals[0])*rdx2SqVolCu[0])*omega+1555848.0*phiC[0]*rdx2SqVolCu[1]+1.1189052e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*phiC[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(624456.0*phiUy[1]-1555848.0*phiC[1])*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiUx[1]-1.1189052e+7*rdx2SqVol[0]*phiC[1]+(847040.3948826761*phiUy[0]+1100685.379244677*phiUx[0]-5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiUx[1]-1.1189052e+7*rdx2SqVolSq[0]*phiC[1]+(38955.5547130316*phiUy[0]+2275410.734360502*phiUx[0]-5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiUx[1]-1555848.0*rdx2SqVolCu[0]*phiC[1]+(366640.5149461798*phiUx[0]-733281.0298923596*bcVals[0])*rdx2SqVolCu[0])*omega+1555848.0*phiC[1]*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]*phiC[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiC[2]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiUx[1]+(2275410.734360502*phiUy[0]+38955.5547130316*phiUx[0]+934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiUx[1]+(1100685.379244677*phiUy[0]+847040.3948826761*phiUx[0]+1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiC[2])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiC[3]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiUx[1]+(90933.0*phiUy[0]+82467.0*phiUx[0]-419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiUx[1]+(82467.0*phiUy[0]+90933.0*phiUx[0]+73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(17043.37994647775*rdx2SqVolCu[1]+381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]+973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1])-2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]-1136585.596333157*rdx2SqVolCu[0])*rho[1]-92496.0*rho[0]*rdx2SqVolCu[1]-2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+((307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(39767.88654178142*rdx2SqVolR4[1]+930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]+2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1]+257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]+747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(115456.0*rdx2SqVolR4[1]+2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]+7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiC[0]-40344.0*phiUy[0])*rdx2SqVolR4[1]+((-310402.557275226*rdx2SqVol[0]*phiUy[1])+10012.98571855568*rdx2SqVol[0]*phiUx[1]+((-945501.0*phiUy[0])-17343.0*phiUx[0]+1170960.0*phiC[0]-416232.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1122732.653974222*rdx2SqVolSq[0]*phiUy[1])+1113691.348758712*rdx2SqVolSq[0]*phiUx[1]+((-2972058.0*phiUy[0])-1286154.0*phiUx[0]+6835740.0*phiC[0]-5155056.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-639329.3979374008*rdx2SqVolCu[0]*phiUy[1])+3757176.73613406*rdx2SqVolCu[0]*phiUx[1]+((-1559817.0*phiUy[0])-4278411.0*phiUx[0]+1.259496e+7*phiC[0]-1.3513464e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1649882.317257809*rdx2SqVolR4[0]*phiUx[1]+((-1873368.0*phiUx[0])+4667544.0*phiC[0]-5588352.0*bcVals[0])*rdx2SqVolR4[0])*omega-40344.0*phiC[0]*rdx2SqVolR4[1]-1170960.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*phiC[0]*rdx2SqVolR4[0]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(85680.0*rdx2SqVol[0]*rdx2SqVolSq[1]+260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-92496.0*rdx2SqVolCu[1])-873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]-825552.0*rdx2SqVolCu[0])*rho[1]-407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0]*rho[0])*omega*volFac+((39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(168555.0*rdx2SqVol[0]*rdx2SqVolCu[1]+606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(20295.0*rdx2SqVol[0]*rdx2SqVolCu[1]+151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1]+1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]+623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiC[1]-40344.0*phiUy[1])*rdx2SqVolR4[1]+((-413649.0*rdx2SqVol[0]*phiUy[1])+109839.0*rdx2SqVol[0]*phiUx[1]+1170960.0*rdx2SqVol[0]*phiC[1]+((-170220.7572154465*phiUy[0])-110142.8429041125*phiUx[0]+560727.200239118*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1048338.0*rdx2SqVolSq[0]*phiUy[1])+1081578.0*rdx2SqVolSq[0]*phiUx[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]+((-615692.1005665087*phiUy[0])-1116705.117163882*phiUx[0]+3464794.435460782*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-464373.0*rdx2SqVolCu[0]*phiUy[1])+2599263.0*rdx2SqVolCu[0]*phiUx[1]+1.259496e+7*rdx2SqVolCu[0]*phiC[1]+((-350599.9924172843*phiUy[0])-2717894.290068507*phiUx[0]+6136988.564971584*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1047816.0*rdx2SqVolR4[0]*phiUx[1]+4667544.0*rdx2SqVolR4[0]*phiC[1]+(2199843.089677078*bcVals[0]-1099921.544838539*phiUx[0])*rdx2SqVolR4[0])*omega-40344.0*phiC[1]*rdx2SqVolR4[1]-1170960.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0]*phiC[1]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = (((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+(17712.0*rdx2SqVolCu[1]+472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0])*rho[2]+((-197904.0*rdx2SqVol[0]*rdx2SqVolSq[1])-600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-21588.28126553848*rho[0]*rdx2SqVolCu[1]-482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiUx[3]+(984.0*rdx2SqVolR4[1]-24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]-659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(3321.0*rdx2SqVol[0]*rdx2SqVolCu[1]+156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]+1873368.0*rdx2SqVolR4[0])*phiUx[2]+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[2]+(45449.01319060734*rdx2SqVolR4[1]+1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]+4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]+3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-72447.0*rdx2SqVol[0]*phiUy[1])+2337.0*rdx2SqVol[0]*phiUx[1]+(52621.43558475006*phiUy[0]-4047.802737288466*phiUx[0]-97147.26569492316*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(287280.0*rdx2SqVolSq[0]*phiUx[1]+(812719.8081306986*phiUy[0]-326193.6644878315*phiUx[0]-973052.2872857346*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(779247.0*rdx2SqVolCu[0]*phiUy[1]+846963.0*rdx2SqVolCu[0]*phiUx[1]+(1901183.83687717*phiUy[0]-946692.2060453439*phiUx[0]-1908983.261663653*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[2])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+(31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1]+270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0])*rho[2]+((-21588.28126553848*rdx2SqVolCu[1])-143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]-204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiUx[3]+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[3]+(39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1]+15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]-376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1]+453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]+1099921.544838539*rdx2SqVolR4[0])*phiUx[2]+(150416.0*rdx2SqVol[0]*rdx2SqVolCu[1]+693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(176754.0528615963*rdx2SqVol[0]*phiUy[1]+25636.08400282695*rdx2SqVol[0]*phiUx[1]+((-39729.0*phiUy[0])-25707.0*phiUx[0]+130872.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(812719.8081306986*rdx2SqVolSq[0]*phiUy[1]+182447.3038660752*rdx2SqVolSq[0]*phiUx[1]+(383040.0*bcVals[0]-191520.0*phiUx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(566001.2949481651*rdx2SqVolCu[0]*phiUy[1]+278113.666120527*rdx2SqVolCu[0]*phiUx[1]+(427329.0*phiUy[0]-304893.0*phiUx[0]-244872.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[3])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-1136585.596333157*rdx2SqVolCu[1])-2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]-742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1]+381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[1]-2773008.0*rho[0]*rdx2SqVolCu[1]-6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(1649882.317257809*rdx2SqVolR4[1]+3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]+1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]+10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1])-1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]-310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+((-5588352.0*rdx2SqVolR4[1])-1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(4667544.0*phiC[0]-1873368.0*phiUy[0])*rdx2SqVolR4[1]+(747388.5837200081*rdx2SqVol[0]*phiUy[1]+1500934.608060923*rdx2SqVol[0]*phiUx[1]+((-4278411.0*phiUy[0])-1559817.0*phiUx[0]+1.259496e+7*phiC[0]+2773008.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(257521.3140693406*rdx2SqVolSq[0]*phiUy[1]+2913967.637637727*rdx2SqVolSq[0]*phiUx[1]+((-1286154.0*phiUy[0])-2972058.0*phiUx[0]+6835740.0*phiC[0]+7782592.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3195.633739964578*rdx2SqVolCu[0]*phiUy[1]+930985.9693223094*rdx2SqVolCu[0]*phiUx[1]+((-17343.0*phiUy[0])-945501.0*phiUx[0]+1170960.0*phiC[0]+2659024.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0]*phiUx[1]+((-40344.0*phiUx[0])+40344.0*phiC[0]+115456.0*bcVals[0])*rdx2SqVolR4[0])*omega-4667544.0*phiC[0]*rdx2SqVolR4[1]-1.259496e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*phiC[0]*rdx2SqVolR4[0]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = (((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-600780.0*rdx2SqVol[0]*rdx2SqVolSq[1])-197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(2773008.0*rdx2SqVolCu[1]+2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]+472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[1]-1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(846963.0*rdx2SqVol[0]*rdx2SqVolCu[1]+287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(779247.0*rdx2SqVol[0]*rdx2SqVolCu[1]-72447.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+((-1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1])-973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]-97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(1873368.0*phiUy[1]-4667544.0*phiC[1])*rdx2SqVolR4[1]+(998973.0*rdx2SqVol[0]*phiUy[1]-1675359.0*rdx2SqVol[0]*phiUx[1]-1.259496e+7*rdx2SqVol[0]*phiC[1]+((-946692.2060453439*phiUy[0])+1901183.83687717*phiUx[0]+3735659.468951633*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(156186.0*rdx2SqVolSq[0]*phiUy[1]-659958.0*rdx2SqVolSq[0]*phiUx[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]+((-326193.6644878315*phiUy[0])+812719.8081306986*phiUx[0]+4193890.830602055*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3321.0*rdx2SqVolCu[0]*phiUy[1]-24363.0*rdx2SqVolCu[0]*phiUx[1]-1170960.0*rdx2SqVolCu[0]*phiC[1]+((-4047.802737288466*phiUy[0])+52621.43558475006*phiUx[0]+1119902.482954654*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+984.0*rdx2SqVolR4[0]*phiUx[1]-40344.0*rdx2SqVolR4[0]*phiC[1]+45449.01319060734*bcVals[0]*rdx2SqVolR4[0])*omega+4667544.0*phiC[1]*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]*phiC[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+((-825552.0*rdx2SqVolCu[1])-2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]-873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0])*rho[2]+(260100.0*rdx2SqVol[0]*rdx2SqVolSq[1]+85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-623288.8754085057*rho[0]*rdx2SqVolCu[1]-1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiUx[3]+(1047816.0*rdx2SqVolR4[1]+2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-464373.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiUx[2]+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[2]+(2199843.089677078*rdx2SqVolR4[1]+6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]+3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]+560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1099921.544838539*phiUy[0]*rdx2SqVolR4[1]+(240705.0*rdx2SqVol[0]*phiUy[1]+337365.0*rdx2SqVol[0]*phiUx[1]+((-2717894.290068507*phiUy[0])-350599.9924172843*phiUx[0]+623288.8754085057*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(151200.0*rdx2SqVolSq[0]*phiUy[1]+606900.0*rdx2SqVolSq[0]*phiUx[1]+((-1116705.117163882*phiUy[0])-615692.1005665087*phiUx[0]+1761980.645523667*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(20295.0*rdx2SqVolCu[0]*phiUy[1]+168555.0*rdx2SqVolCu[0]*phiUx[1]+((-110142.8429041125*phiUy[0])-170220.7572154465*phiUx[0]+522469.6620015367*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[2]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+((-204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1])-143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0])*rho[2]+(623288.8754085057*rdx2SqVolCu[1]+270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]+31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiUx[3]+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[3]+(278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1]+182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]+25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1]+812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]+176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+((-244872.0*rdx2SqVol[0]*rdx2SqVolCu[1])+383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1099921.544838539*phiUy[1]*rdx2SqVolR4[1]+(1661713.956324312*rdx2SqVol[0]*phiUy[1]-376570.3622259722*rdx2SqVol[0]*phiUx[1]+((-304893.0*phiUy[0])+427329.0*phiUx[0]+839664.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(453260.3758326994*rdx2SqVolSq[0]*phiUy[1]+15016.88050162216*rdx2SqVolSq[0]*phiUx[1]+(693600.0*bcVals[0]-191520.0*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(21091.18268376621*rdx2SqVolCu[0]*phiUy[1]+39956.68007980643*rdx2SqVolCu[0]*phiUx[1]+((-25707.0*phiUy[0])-39729.0*phiUx[0]+150416.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[3])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiUy[0]-40344.0*phiC[0])*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiUx[1]+(293355.0*phiUy[0]+52029.0*phiUx[0]-345384.0*phiC[0]-92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+(52029.0*phiUy[0]+293355.0*phiUx[0]-345384.0*phiC[0]-828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+(40344.0*phiUx[0]-40344.0*phiC[0]-115456.0*bcVals[0])*rdx2SqVolCu[0])*omega+40344.0*phiC[0]*rdx2SqVolCu[1]+345384.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*phiC[0]*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(121032.0*phiC[1]-121032.0*phiUy[1])*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiUx[1]+1036152.0*rdx2SqVol[0]*phiC[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiUx[0]-373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiUx[1]+1036152.0*rdx2SqVolSq[0]*phiC[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiUx[0]-1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiUx[1]+121032.0*rdx2SqVolCu[0]*phiC[1]-136347.039571822*bcVals[0]*rdx2SqVolCu[0])*omega-121032.0*phiC[1]*rdx2SqVolCu[1]-1036152.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0]*phiC[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[2]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiUx[1]+((-36430.22463559618*phiUy[0])+36430.22463559618*phiUx[0]-64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiUx[1]+((-190246.7286525579*phiUy[0])+190246.7286525579*phiUx[0]-583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[2]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[3]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUy[2]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiUx[1]+(44403.0*phiUy[0]-44403.0*phiUx[0]-87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiUx[1]+((-44403.0*phiUy[0])+44403.0*phiUx[0]-168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(1555848.0*phiC[0]-624456.0*phiLy[0])*rdx2SqVolCu[1]+(100062.3072040616*rdx2SqVol[0]*phiUx[1]-1544603.07302135*rdx2SqVol[0]*phiLy[1]+((-173313.0*phiUx[0])-3368931.0*phiLy[0]+1.1189052e+7*phiC[0]-4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2951378.203030407*rdx2SqVolSq[0]*phiUx[1]-71036.59977082233*rdx2SqVolSq[0]*phiLy[1]+((-3368931.0*phiUx[0])-173313.0*phiLy[0]+1.1189052e+7*phiC[0]-1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+((-624456.0*phiUx[0])+1555848.0*phiC[0]-1862784.0*bcVals[0])*rdx2SqVolCu[0])*omega-1555848.0*phiC[0]*rdx2SqVolCu[1]-1.1189052e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*phiC[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(1555848.0*phiC[1]-624456.0*phiLy[1])*rdx2SqVolCu[1]+(1097649.0*rdx2SqVol[0]*phiUx[1]-722367.0*rdx2SqVol[0]*phiLy[1]+1.1189052e+7*rdx2SqVol[0]*phiC[1]+((-1100685.379244677*phiUx[0])-847040.3948826761*phiLy[0]+5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2182239.0*rdx2SqVolSq[0]*phiUx[1]-51597.0*rdx2SqVolSq[0]*phiLy[1]+1.1189052e+7*rdx2SqVolSq[0]*phiC[1]+((-2275410.734360502*phiUx[0])-38955.5547130316*phiLy[0]+5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiUx[1]+1555848.0*rdx2SqVolCu[0]*phiC[1]+(733281.0298923596*bcVals[0]-366640.5149461798*phiUx[0])*rdx2SqVolCu[0])*omega-1555848.0*phiC[1]*rdx2SqVolCu[1]-1.1189052e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0]*phiC[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiC[2]-366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+(22491.0*rdx2SqVol[0]*phiUx[1]-497457.0*rdx2SqVol[0]*phiLy[1]+((-38955.5547130316*phiUx[0])-2275410.734360502*phiLy[0]-934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(757809.0*rdx2SqVolSq[0]*phiUx[1]-451143.0*rdx2SqVolSq[0]*phiLy[1]+((-847040.3948826761*phiUx[0])-1100685.379244677*phiLy[0]-1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiC[2])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiC[3]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[3]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+(82239.50439417786*rdx2SqVol[0]*phiUx[1]-474351.5585164657*rdx2SqVol[0]*phiLy[1]+((-82467.0*phiUx[0])-90933.0*phiLy[0]+419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(82946.18112366594*rdx2SqVolSq[0]*phiUx[1]-109228.3200777161*rdx2SqVolSq[0]*phiLy[1]+((-90933.0*phiUx[0])-82467.0*phiLy[0]-73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(17043.37994647775*rdx2SqVolCu[1]+381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]+973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1]+2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[1]+92496.0*rho[0]*rdx2SqVolCu[1]+2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(115456.0*rdx2SqVolR4[1]+2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]+7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1]+257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]+747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(39767.88654178142*rdx2SqVolR4[1]+930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]+2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(40344.0*phiLy[0]-40344.0*phiC[0])*rdx2SqVolR4[1]+((-10012.98571855568*rdx2SqVol[0]*phiUx[1])+310402.557275226*rdx2SqVol[0]*phiLy[1]+(17343.0*phiUx[0]+945501.0*phiLy[0]-1170960.0*phiC[0]+416232.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1113691.348758712*rdx2SqVolSq[0]*phiUx[1])+1122732.653974222*rdx2SqVolSq[0]*phiLy[1]+(1286154.0*phiUx[0]+2972058.0*phiLy[0]-6835740.0*phiC[0]+5155056.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-3757176.73613406*rdx2SqVolCu[0]*phiUx[1])+639329.3979374008*rdx2SqVolCu[0]*phiLy[1]+(4278411.0*phiUx[0]+1559817.0*phiLy[0]-1.259496e+7*phiC[0]+1.3513464e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0]*phiUx[1]+(1873368.0*phiUx[0]-4667544.0*phiC[0]+5588352.0*bcVals[0])*rdx2SqVolR4[0])*omega+40344.0*phiC[0]*rdx2SqVolR4[1]+1170960.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*phiC[0]*rdx2SqVolR4[0])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = (((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(85680.0*rdx2SqVol[0]*rdx2SqVolSq[1]+260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(92496.0*rdx2SqVolCu[1]+873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[1]+407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+(39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1]+1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]+623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(20295.0*rdx2SqVol[0]*rdx2SqVolCu[1]+151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(168555.0*rdx2SqVol[0]*rdx2SqVolCu[1]+606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(40344.0*phiLy[1]-40344.0*phiC[1])*rdx2SqVolR4[1]+((-109839.0*rdx2SqVol[0]*phiUx[1])+413649.0*rdx2SqVol[0]*phiLy[1]-1170960.0*rdx2SqVol[0]*phiC[1]+(110142.8429041125*phiUx[0]+170220.7572154465*phiLy[0]-560727.200239118*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1081578.0*rdx2SqVolSq[0]*phiUx[1])+1048338.0*rdx2SqVolSq[0]*phiLy[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]+(1116705.117163882*phiUx[0]+615692.1005665087*phiLy[0]-3464794.435460782*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2599263.0*rdx2SqVolCu[0]*phiUx[1])+464373.0*rdx2SqVolCu[0]*phiLy[1]-1.259496e+7*rdx2SqVolCu[0]*phiC[1]+(2717894.290068507*phiUx[0]+350599.9924172843*phiLy[0]-6136988.564971584*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0]*phiUx[1]-4667544.0*rdx2SqVolR4[0]*phiC[1]+(1099921.544838539*phiUx[0]-2199843.089677078*bcVals[0])*rdx2SqVolR4[0])*omega+40344.0*phiC[1]*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]*phiC[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = (((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+(17712.0*rdx2SqVolCu[1]+472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0])*rho[2]+(197904.0*rdx2SqVol[0]*rdx2SqVolSq[1]+600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+21588.28126553848*rho[0]*rdx2SqVolCu[1]+482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiUx[3]+(72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(45449.01319060734*rdx2SqVolR4[1]+1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]+4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]+3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(3321.0*rdx2SqVol[0]*rdx2SqVolCu[1]+156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]+1873368.0*rdx2SqVolR4[0])*phiUx[2]+(984.0*rdx2SqVolR4[1]-24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]-659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[2]+((-2337.0*rdx2SqVol[0]*phiUx[1])+72447.0*rdx2SqVol[0]*phiLy[1]+(4047.802737288466*phiUx[0]-52621.43558475006*phiLy[0]+97147.26569492316*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((326193.6644878315*phiUx[0]-812719.8081306986*phiLy[0]+973052.2872857346*bcVals[0])*rdx2SqVolSq[0]-287280.0*rdx2SqVolSq[0]*phiUx[1])*rdx2SqVolSq[1]+((-846963.0*rdx2SqVolCu[0]*phiUx[1])-779247.0*rdx2SqVolCu[0]*phiLy[1]+(946692.2060453439*phiUx[0]-1901183.83687717*phiLy[0]+1908983.261663653*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[2])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+(31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1]+270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0])*rho[2]+(21588.28126553848*rdx2SqVolCu[1]+143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]+204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiUx[3]+(984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[3]+(150416.0*rdx2SqVol[0]*rdx2SqVolCu[1]+693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1]+453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]+1099921.544838539*rdx2SqVolR4[0])*phiUx[2]+(39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1]+15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]-376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-25636.08400282695*rdx2SqVol[0]*phiUx[1])-176754.0528615963*rdx2SqVol[0]*phiLy[1]+(25707.0*phiUx[0]+39729.0*phiLy[0]-130872.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-182447.3038660752*rdx2SqVolSq[0]*phiUx[1])-812719.8081306986*rdx2SqVolSq[0]*phiLy[1]+(191520.0*phiUx[0]-383040.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-278113.666120527*rdx2SqVolCu[0]*phiUx[1])-566001.2949481651*rdx2SqVolCu[0]*phiLy[1]+(304893.0*phiUx[0]-427329.0*phiLy[0]+244872.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[3])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-1136585.596333157*rdx2SqVolCu[1])-2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]-742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1])-381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]-17043.37994647775*rdx2SqVolCu[0])*rho[1]+2773008.0*rho[0]*rdx2SqVolCu[1]+6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+((615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(5588352.0*rdx2SqVolR4[1]+1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1])-1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]-310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(1649882.317257809*rdx2SqVolR4[1]+3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]+1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]+10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(1873368.0*phiLy[0]-4667544.0*phiC[0])*rdx2SqVolR4[1]+((-1500934.608060923*rdx2SqVol[0]*phiUx[1])-747388.5837200081*rdx2SqVol[0]*phiLy[1]+(1559817.0*phiUx[0]+4278411.0*phiLy[0]-1.259496e+7*phiC[0]-2773008.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-2913967.637637727*rdx2SqVolSq[0]*phiUx[1])-257521.3140693406*rdx2SqVolSq[0]*phiLy[1]+(2972058.0*phiUx[0]+1286154.0*phiLy[0]-6835740.0*phiC[0]-7782592.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-930985.9693223094*rdx2SqVolCu[0]*phiUx[1])-3195.633739964578*rdx2SqVolCu[0]*phiLy[1]+(945501.0*phiUx[0]+17343.0*phiLy[0]-1170960.0*phiC[0]-2659024.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolR4[0]*phiUx[1]+(40344.0*phiUx[0]-40344.0*phiC[0]-115456.0*bcVals[0])*rdx2SqVolR4[0])*omega+4667544.0*phiC[0]*rdx2SqVolR4[1]+1.259496e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*phiC[0]*rdx2SqVolR4[0])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-600780.0*rdx2SqVol[0]*rdx2SqVolSq[1])-197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-2773008.0*rdx2SqVolCu[1])-2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]-472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]-17712.0*rdx2SqVolCu[0])*rho[1]+1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1]+973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]+97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(779247.0*rdx2SqVol[0]*rdx2SqVolCu[1]-72447.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(846963.0*rdx2SqVol[0]*rdx2SqVolCu[1]+287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(4667544.0*phiC[1]-1873368.0*phiLy[1])*rdx2SqVolR4[1]+(1675359.0*rdx2SqVol[0]*phiUx[1]-998973.0*rdx2SqVol[0]*phiLy[1]+1.259496e+7*rdx2SqVol[0]*phiC[1]+((-1901183.83687717*phiUx[0])+946692.2060453439*phiLy[0]-3735659.468951633*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(659958.0*rdx2SqVolSq[0]*phiUx[1]-156186.0*rdx2SqVolSq[0]*phiLy[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]+((-812719.8081306986*phiUx[0])+326193.6644878315*phiLy[0]-4193890.830602055*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(24363.0*rdx2SqVolCu[0]*phiUx[1]-3321.0*rdx2SqVolCu[0]*phiLy[1]+1170960.0*rdx2SqVolCu[0]*phiC[1]+((-52621.43558475006*phiUx[0])+4047.802737288466*phiLy[0]-1119902.482954654*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-984.0*rdx2SqVolR4[0]*phiUx[1]+40344.0*rdx2SqVolR4[0]*phiC[1]-45449.01319060734*bcVals[0]*rdx2SqVolR4[0])*omega-4667544.0*phiC[1]*rdx2SqVolR4[1]-1.259496e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0]*phiC[1]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+((-825552.0*rdx2SqVolCu[1])-2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]-873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0])*rho[2]+((-260100.0*rdx2SqVol[0]*rdx2SqVolSq[1])-85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+623288.8754085057*rho[0]*rdx2SqVolCu[1]+1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiUx[3]+((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-2199843.089677078*rdx2SqVolR4[1])-6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]-3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]-560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-464373.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiUx[2]+(1047816.0*rdx2SqVolR4[1]+2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[2]+1099921.544838539*phiLy[0]*rdx2SqVolR4[1]+((-337365.0*rdx2SqVol[0]*phiUx[1])-240705.0*rdx2SqVol[0]*phiLy[1]+(350599.9924172843*phiUx[0]+2717894.290068507*phiLy[0]-623288.8754085057*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-606900.0*rdx2SqVolSq[0]*phiUx[1])-151200.0*rdx2SqVolSq[0]*phiLy[1]+(615692.1005665087*phiUx[0]+1116705.117163882*phiLy[0]-1761980.645523667*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-168555.0*rdx2SqVolCu[0]*phiUx[1])-20295.0*rdx2SqVolCu[0]*phiLy[1]+(170220.7572154465*phiUx[0]+110142.8429041125*phiLy[0]-522469.6620015367*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[2]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+((-204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1])-143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0])*rho[2]+((-623288.8754085057*rdx2SqVolCu[1])-270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]-31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiUx[3]+((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[3]+(244872.0*rdx2SqVol[0]*rdx2SqVolCu[1]-383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1]+812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]+176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1]+182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]+25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-1099921.544838539*phiLy[1]*rdx2SqVolR4[1]+(376570.3622259722*rdx2SqVol[0]*phiUx[1]-1661713.956324312*rdx2SqVol[0]*phiLy[1]+((-427329.0*phiUx[0])+304893.0*phiLy[0]-839664.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-15016.88050162216*rdx2SqVolSq[0]*phiUx[1])-453260.3758326994*rdx2SqVolSq[0]*phiLy[1]+(191520.0*phiLy[0]-693600.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-39956.68007980643*rdx2SqVolCu[0]*phiUx[1])-21091.18268376621*rdx2SqVolCu[0]*phiLy[1]+(39729.0*phiUx[0]+25707.0*phiLy[0]-150416.0*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[3])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(40344.0*phiC[0]-40344.0*phiLy[0])*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiUx[1]+50064.92859277839*rdx2SqVol[0]*phiLy[1]+((-52029.0*phiUx[0])-293355.0*phiLy[0]+345384.0*phiC[0]+92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+((-293355.0*phiUx[0])-52029.0*phiLy[0]+345384.0*phiC[0]+828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+((-40344.0*phiUx[0])+40344.0*phiC[0]+115456.0*bcVals[0])*rdx2SqVolCu[0])*omega-40344.0*phiC[0]*rdx2SqVolCu[1]-345384.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-345384.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-40344.0*phiC[0]*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(121032.0*phiLy[1]-121032.0*phiC[1])*rdx2SqVolCu[1]+((-167649.0*rdx2SqVol[0]*phiUx[1])+221031.0*rdx2SqVol[0]*phiLy[1]-1036152.0*rdx2SqVol[0]*phiC[1]+(190246.7286525579*phiUx[0]-190246.7286525579*phiLy[0]+373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-11367.0*rdx2SqVolSq[0]*phiUx[1])+29889.0*rdx2SqVolSq[0]*phiLy[1]-1036152.0*rdx2SqVolSq[0]*phiC[1]+(36430.22463559618*phiUx[0]-36430.22463559618*phiLy[0]+1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiUx[1]-121032.0*rdx2SqVolCu[0]*phiC[1]+136347.039571822*bcVals[0]*rdx2SqVolCu[0])*omega+121032.0*phiC[1]*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]*phiC[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[2]+(35055.0*rdx2SqVol[0]*phiUx[1]+35055.0*rdx2SqVol[0]*phiLy[1]+((-36430.22463559618*phiUx[0])+36430.22463559618*phiLy[0]+64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(188385.0*rdx2SqVolSq[0]*phiUx[1]-35055.0*rdx2SqVolSq[0]*phiLy[1]+((-190246.7286525579*phiUx[0])+190246.7286525579*phiLy[0]+583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[2]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[3]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLy[2]+((-39128.75979378851*rdx2SqVol[0]*phiUx[1])-190246.7286525579*rdx2SqVol[0]*phiLy[1]+(44403.0*phiUx[0]-44403.0*phiLy[0]+87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(44657.46597154836*rdx2SqVolSq[0]*phiUx[1]-36430.22463559618*rdx2SqVolSq[0]*phiLy[1]+((-44403.0*phiUx[0])+44403.0*phiLy[0]+168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(1555848.0*phiC[0]-624456.0*phiUy[0])*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiLx[1]-4159512.0*rdx2SqVol[0]*bcVals[1]+((-3368931.0*phiUy[0])-173313.0*phiLx[0]+1.1189052e+7*phiC[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiLx[1]-1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+((-173313.0*phiUy[0])-3368931.0*phiLx[0]+1.1189052e+7*phiC[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiLx[1]-1862784.0*rdx2SqVolCu[0]*bcVals[1]+(1555848.0*phiC[0]-624456.0*phiLx[0])*rdx2SqVolCu[0])*omega-1555848.0*phiC[0]*rdx2SqVolCu[1]-1.1189052e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*phiC[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(624456.0*phiUy[1]-1555848.0*phiC[1])*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiLx[1]-1.1189052e+7*rdx2SqVol[0]*phiC[1]+5603489.203427449*rdx2SqVol[0]*bcVals[1]+((-847040.3948826761*phiUy[0])-1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiLx[1]-1.1189052e+7*rdx2SqVolSq[0]*phiC[1]+5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+((-38955.5547130316*phiUy[0])-2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiLx[1]-1555848.0*rdx2SqVolCu[0]*phiC[1]+733281.0298923596*rdx2SqVolCu[0]*bcVals[1]-366640.5149461798*phiLx[0]*rdx2SqVolCu[0])*omega+1555848.0*phiC[1]*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]*phiC[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiC[2]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiLx[1]-934933.3131127583*rdx2SqVol[0]*bcVals[1]+((-2275410.734360502*phiUy[0])-38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiLx[1]-1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+((-1100685.379244677*phiUy[0])-847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiC[2]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiC[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiLx[1]+419832.0*rdx2SqVol[0]*bcVals[1]+((-90933.0*phiUy[0])-82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiLx[1]-73032.0*rdx2SqVolSq[0]*bcVals[1]+((-82467.0*phiUy[0])-90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-17043.37994647775*rdx2SqVolCu[1])-381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]-973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1])-2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]-1136585.596333157*rdx2SqVolCu[0])*rho[1]+92496.0*rho[0]*rdx2SqVolCu[1]+2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+2773008.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+((307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-39767.88654178142*rdx2SqVolR4[1])-930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]-2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1])-257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]-747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-115456.0*rdx2SqVolR4[1])-2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]-7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiUy[0]-40344.0*phiC[0])*rdx2SqVolR4[1]+((-310402.557275226*rdx2SqVol[0]*phiUy[1])+10012.98571855568*rdx2SqVol[0]*phiLx[1]+416232.0*rdx2SqVol[0]*bcVals[1]+(945501.0*phiUy[0]+17343.0*phiLx[0]-1170960.0*phiC[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1122732.653974222*rdx2SqVolSq[0]*phiUy[1])+1113691.348758712*rdx2SqVolSq[0]*phiLx[1]+5155056.0*rdx2SqVolSq[0]*bcVals[1]+(2972058.0*phiUy[0]+1286154.0*phiLx[0]-6835740.0*phiC[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-639329.3979374008*rdx2SqVolCu[0]*phiUy[1])+3757176.73613406*rdx2SqVolCu[0]*phiLx[1]+1.3513464e+7*rdx2SqVolCu[0]*bcVals[1]+(1559817.0*phiUy[0]+4278411.0*phiLx[0]-1.259496e+7*phiC[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1649882.317257809*rdx2SqVolR4[0]*phiLx[1]+5588352.0*rdx2SqVolR4[0]*bcVals[1]+(1873368.0*phiLx[0]-4667544.0*phiC[0])*rdx2SqVolR4[0])*omega+40344.0*phiC[0]*rdx2SqVolR4[1]+1170960.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*phiC[0]*rdx2SqVolR4[0])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-85680.0*rdx2SqVol[0]*rdx2SqVolSq[1])-260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-92496.0*rdx2SqVolCu[1])-873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]-825552.0*rdx2SqVolCu[0])*rho[1]+407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+623288.8754085057*rdx2SqVolCu[0]*rho[0])*omega*volFac+((39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-168555.0*rdx2SqVol[0]*rdx2SqVolCu[1])-606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-20295.0*rdx2SqVol[0]*rdx2SqVolCu[1])-151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1])-1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]-623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiC[1]-40344.0*phiUy[1])*rdx2SqVolR4[1]+((-413649.0*rdx2SqVol[0]*phiUy[1])+109839.0*rdx2SqVol[0]*phiLx[1]+1170960.0*rdx2SqVol[0]*phiC[1]-560727.200239118*rdx2SqVol[0]*bcVals[1]+(170220.7572154465*phiUy[0]+110142.8429041125*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-1048338.0*rdx2SqVolSq[0]*phiUy[1])+1081578.0*rdx2SqVolSq[0]*phiLx[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]-3464794.435460782*rdx2SqVolSq[0]*bcVals[1]+(615692.1005665087*phiUy[0]+1116705.117163882*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-464373.0*rdx2SqVolCu[0]*phiUy[1])+2599263.0*rdx2SqVolCu[0]*phiLx[1]+1.259496e+7*rdx2SqVolCu[0]*phiC[1]-6136988.564971584*rdx2SqVolCu[0]*bcVals[1]+(350599.9924172843*phiUy[0]+2717894.290068507*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1047816.0*rdx2SqVolR4[0]*phiLx[1]+4667544.0*rdx2SqVolR4[0]*phiC[1]-2199843.089677078*rdx2SqVolR4[0]*bcVals[1]+1099921.544838539*phiLx[0]*rdx2SqVolR4[0])*omega-40344.0*phiC[1]*rdx2SqVolR4[1]-1170960.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0]*phiC[1]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+((-17712.0*rdx2SqVolCu[1])-472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0])*rho[2]+((-197904.0*rdx2SqVol[0]*rdx2SqVolSq[1])-600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+21588.28126553848*rho[0]*rdx2SqVolCu[1]+482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiLx[3]+((-984.0*rdx2SqVolR4[1])+24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]+659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-3321.0*rdx2SqVol[0]*rdx2SqVolCu[1])-156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1873368.0*rdx2SqVolR4[0])*phiLx[2]+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[2]+((-45449.01319060734*rdx2SqVolR4[1])-1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]-4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-72447.0*rdx2SqVol[0]*phiUy[1])+2337.0*rdx2SqVol[0]*phiLx[1]+97147.26569492316*rdx2SqVol[0]*bcVals[1]+(4047.802737288466*phiLx[0]-52621.43558475006*phiUy[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(287280.0*rdx2SqVolSq[0]*phiLx[1]+973052.2872857346*rdx2SqVolSq[0]*bcVals[1]+(326193.6644878315*phiLx[0]-812719.8081306986*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(779247.0*rdx2SqVolCu[0]*phiUy[1]+846963.0*rdx2SqVolCu[0]*phiLx[1]+1908983.261663653*rdx2SqVolCu[0]*bcVals[1]+(946692.2060453439*phiLx[0]-1901183.83687717*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[2]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+((-31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1])-270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0])*rho[2]+((-21588.28126553848*rdx2SqVolCu[1])-143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]-204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiLx[3]+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[3]+((-39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1])-15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]+376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1])-453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]-1099921.544838539*rdx2SqVolR4[0])*phiLx[2]+((-150416.0*rdx2SqVol[0]*rdx2SqVolCu[1])-693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(176754.0528615963*rdx2SqVol[0]*phiUy[1]+25636.08400282695*rdx2SqVol[0]*phiLx[1]-130872.0*rdx2SqVol[0]*bcVals[1]+(39729.0*phiUy[0]+25707.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(812719.8081306986*rdx2SqVolSq[0]*phiUy[1]+182447.3038660752*rdx2SqVolSq[0]*phiLx[1]-383040.0*rdx2SqVolSq[0]*bcVals[1]+191520.0*phiLx[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(566001.2949481651*rdx2SqVolCu[0]*phiUy[1]+278113.666120527*rdx2SqVolCu[0]*phiLx[1]+244872.0*rdx2SqVolCu[0]*bcVals[1]+(304893.0*phiLx[0]-427329.0*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[3])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1136585.596333157*rdx2SqVolCu[1]+2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]+742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1]+381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[1]+2773008.0*rho[0]*rdx2SqVolCu[1]+6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-1649882.317257809*rdx2SqVolR4[1])-3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]-1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]-10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1]+1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]+310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(5588352.0*rdx2SqVolR4[1]+1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(1873368.0*phiUy[0]-4667544.0*phiC[0])*rdx2SqVolR4[1]+(747388.5837200081*rdx2SqVol[0]*phiUy[1]+1500934.608060923*rdx2SqVol[0]*phiLx[1]+2773008.0*rdx2SqVol[0]*bcVals[1]+(4278411.0*phiUy[0]+1559817.0*phiLx[0]-1.259496e+7*phiC[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(257521.3140693406*rdx2SqVolSq[0]*phiUy[1]+2913967.637637727*rdx2SqVolSq[0]*phiLx[1]+7782592.0*rdx2SqVolSq[0]*bcVals[1]+(1286154.0*phiUy[0]+2972058.0*phiLx[0]-6835740.0*phiC[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3195.633739964578*rdx2SqVolCu[0]*phiUy[1]+930985.9693223094*rdx2SqVolCu[0]*phiLx[1]+2659024.0*rdx2SqVolCu[0]*bcVals[1]+(17343.0*phiUy[0]+945501.0*phiLx[0]-1170960.0*phiC[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0]*phiLx[1]+115456.0*rdx2SqVolR4[0]*bcVals[1]+(40344.0*phiLx[0]-40344.0*phiC[0])*rdx2SqVolR4[0])*omega+4667544.0*phiC[0]*rdx2SqVolR4[1]+1.259496e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*phiC[0]*rdx2SqVolR4[0])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = (((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(600780.0*rdx2SqVol[0]*rdx2SqVolSq[1]+197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(2773008.0*rdx2SqVolCu[1]+2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]+472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[1]+1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-846963.0*rdx2SqVol[0]*rdx2SqVolCu[1])-287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(72447.0*rdx2SqVolCu[0]*rdx2SqVol[1]-779247.0*rdx2SqVol[0]*rdx2SqVolCu[1])*phiLx[2]+(1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1]+973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]+97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+(1873368.0*phiUy[1]-4667544.0*phiC[1])*rdx2SqVolR4[1]+(998973.0*rdx2SqVol[0]*phiUy[1]-1675359.0*rdx2SqVol[0]*phiLx[1]-1.259496e+7*rdx2SqVol[0]*phiC[1]+3735659.468951633*rdx2SqVol[0]*bcVals[1]+(946692.2060453439*phiUy[0]-1901183.83687717*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(156186.0*rdx2SqVolSq[0]*phiUy[1]-659958.0*rdx2SqVolSq[0]*phiLx[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]+4193890.830602055*rdx2SqVolSq[0]*bcVals[1]+(326193.6644878315*phiUy[0]-812719.8081306986*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(3321.0*rdx2SqVolCu[0]*phiUy[1]-24363.0*rdx2SqVolCu[0]*phiLx[1]-1170960.0*rdx2SqVolCu[0]*phiC[1]+1119902.482954654*rdx2SqVolCu[0]*bcVals[1]+(4047.802737288466*phiUy[0]-52621.43558475006*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+984.0*rdx2SqVolR4[0]*phiLx[1]-40344.0*rdx2SqVolR4[0]*phiC[1]+45449.01319060734*rdx2SqVolR4[0]*bcVals[1])*omega+4667544.0*phiC[1]*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]*phiC[1])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = (((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+(825552.0*rdx2SqVolCu[1]+2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]+873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0])*rho[2]+(260100.0*rdx2SqVol[0]*rdx2SqVolSq[1]+85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+623288.8754085057*rho[0]*rdx2SqVolCu[1]+1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiLx[3]+((-1047816.0*rdx2SqVolR4[1])-2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]-1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(464373.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiLx[2]+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[2]+((-2199843.089677078*rdx2SqVolR4[1])-6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]-3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]-560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1099921.544838539*phiUy[0]*rdx2SqVolR4[1]+(240705.0*rdx2SqVol[0]*phiUy[1]+337365.0*rdx2SqVol[0]*phiLx[1]+623288.8754085057*rdx2SqVol[0]*bcVals[1]+(2717894.290068507*phiUy[0]+350599.9924172843*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(151200.0*rdx2SqVolSq[0]*phiUy[1]+606900.0*rdx2SqVolSq[0]*phiLx[1]+1761980.645523667*rdx2SqVolSq[0]*bcVals[1]+(1116705.117163882*phiUy[0]+615692.1005665087*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(20295.0*rdx2SqVolCu[0]*phiUy[1]+168555.0*rdx2SqVolCu[0]*phiLx[1]+522469.6620015367*rdx2SqVolCu[0]*bcVals[1]+(110142.8429041125*phiUy[0]+170220.7572154465*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[2])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+(204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1]+143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0])*rho[2]+(623288.8754085057*rdx2SqVolCu[1]+270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]+31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiLx[3]+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[3]+((-278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1])-182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]-25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1])-812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]-176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(244872.0*rdx2SqVol[0]*rdx2SqVolCu[1]-383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+1099921.544838539*phiUy[1]*rdx2SqVolR4[1]+(1661713.956324312*rdx2SqVol[0]*phiUy[1]-376570.3622259722*rdx2SqVol[0]*phiLx[1]+839664.0*rdx2SqVol[0]*bcVals[1]+(304893.0*phiUy[0]-427329.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(453260.3758326994*rdx2SqVolSq[0]*phiUy[1]+15016.88050162216*rdx2SqVolSq[0]*phiLx[1]+693600.0*rdx2SqVolSq[0]*bcVals[1]+191520.0*phiUy[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(21091.18268376621*rdx2SqVolCu[0]*phiUy[1]+39956.68007980643*rdx2SqVolCu[0]*phiLx[1]+150416.0*rdx2SqVolCu[0]*bcVals[1]+(25707.0*phiUy[0]+39729.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[3])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiC[0]-40344.0*phiUy[0])*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiLx[1]-92496.0*rdx2SqVol[0]*bcVals[1]+((-293355.0*phiUy[0])-52029.0*phiLx[0]+345384.0*phiC[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiLx[1]-828720.0*rdx2SqVolSq[0]*bcVals[1]+((-52029.0*phiUy[0])-293355.0*phiLx[0]+345384.0*phiC[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiLx[1]-115456.0*rdx2SqVolCu[0]*bcVals[1]+(40344.0*phiC[0]-40344.0*phiLx[0])*rdx2SqVolCu[0])*omega-40344.0*phiC[0]*rdx2SqVolCu[1]-345384.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-345384.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-40344.0*phiC[0]*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(121032.0*phiC[1]-121032.0*phiUy[1])*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiLx[1]+1036152.0*rdx2SqVol[0]*phiC[1]-373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiLx[1]+1036152.0*rdx2SqVolSq[0]*phiC[1]-1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiLx[1]+121032.0*rdx2SqVolCu[0]*phiC[1]-136347.039571822*rdx2SqVolCu[0]*bcVals[1])*omega-121032.0*phiC[1]*rdx2SqVolCu[1]-1036152.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0]*phiC[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[2]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiLx[1]-64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiLx[1]-583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[2])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiLx[1]-87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiLx[1]-168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiUy[0]-44403.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(624456.0*phiLy[0]-1555848.0*phiC[0])*rdx2SqVolCu[1]+((-1544603.07302135*rdx2SqVol[0]*phiLy[1])+100062.3072040616*rdx2SqVol[0]*phiLx[1]+4159512.0*rdx2SqVol[0]*bcVals[1]+(3368931.0*phiLy[0]+173313.0*phiLx[0]-1.1189052e+7*phiC[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-71036.59977082233*rdx2SqVolSq[0]*phiLy[1])+2951378.203030407*rdx2SqVolSq[0]*phiLx[1]+1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+(173313.0*phiLy[0]+3368931.0*phiLx[0]-1.1189052e+7*phiC[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiLx[1]+1862784.0*rdx2SqVolCu[0]*bcVals[1]+(624456.0*phiLx[0]-1555848.0*phiC[0])*rdx2SqVolCu[0])*omega+1555848.0*phiC[0]*rdx2SqVolCu[1]+1.1189052e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*phiC[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(1555848.0*phiC[1]-624456.0*phiLy[1])*rdx2SqVolCu[1]+((-722367.0*rdx2SqVol[0]*phiLy[1])+1097649.0*rdx2SqVol[0]*phiLx[1]+1.1189052e+7*rdx2SqVol[0]*phiC[1]-5603489.203427449*rdx2SqVol[0]*bcVals[1]+(847040.3948826761*phiLy[0]+1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-51597.0*rdx2SqVolSq[0]*phiLy[1])+2182239.0*rdx2SqVolSq[0]*phiLx[1]+1.1189052e+7*rdx2SqVolSq[0]*phiC[1]-5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+(38955.5547130316*phiLy[0]+2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiLx[1]+1555848.0*rdx2SqVolCu[0]*phiC[1]-733281.0298923596*rdx2SqVolCu[0]*bcVals[1]+366640.5149461798*phiLx[0]*rdx2SqVolCu[0])*omega-1555848.0*phiC[1]*rdx2SqVolCu[1]-1.1189052e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0]*phiC[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiC[2]+366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+((-497457.0*rdx2SqVol[0]*phiLy[1])+22491.0*rdx2SqVol[0]*phiLx[1]+934933.3131127583*rdx2SqVol[0]*bcVals[1]+(2275410.734360502*phiLy[0]+38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-451143.0*rdx2SqVolSq[0]*phiLy[1])+757809.0*rdx2SqVolSq[0]*phiLx[1]+1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+(1100685.379244677*phiLy[0]+847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiC[2]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiC[3]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+((-474351.5585164657*rdx2SqVol[0]*phiLy[1])+82239.50439417786*rdx2SqVol[0]*phiLx[1]-419832.0*rdx2SqVol[0]*bcVals[1]+(90933.0*phiLy[0]+82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-109228.3200777161*rdx2SqVolSq[0]*phiLy[1])+82946.18112366594*rdx2SqVolSq[0]*phiLx[1]+73032.0*rdx2SqVolSq[0]*bcVals[1]+(82467.0*phiLy[0]+90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((156240.0*rdx2SqVol[0]*rdx2SqVolSq[1]+474300.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-17043.37994647775*rdx2SqVolCu[1])-381189.7417297584*rdx2SqVol[0]*rdx2SqVolSq[1]-973862.8870636768*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(742259.9812787967*rdx2SqVol[0]*rdx2SqVolSq[1]+2574069.987160411*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[1]-92496.0*rho[0]*rdx2SqVolCu[1]-2145504.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-6470652.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+((307365.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+615195.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-1845.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-668655.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-115456.0*rdx2SqVolR4[1])-2659024.0*rdx2SqVol[0]*rdx2SqVolCu[1]-7782592.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2773008.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-39767.88654178142*rdx2SqVolR4[1])-930985.9693223094*rdx2SqVol[0]*rdx2SqVolCu[1]-2913967.637637727*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1500934.608060923*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-3195.633739964578*rdx2SqVol[0]*rdx2SqVolCu[1])-257521.3140693406*rdx2SqVolSq[0]*rdx2SqVolSq[1]-747388.5837200081*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(40344.0*phiC[0]-40344.0*phiLy[0])*rdx2SqVolR4[1]+(310402.557275226*rdx2SqVol[0]*phiLy[1]-10012.98571855568*rdx2SqVol[0]*phiLx[1]-416232.0*rdx2SqVol[0]*bcVals[1]+((-945501.0*phiLy[0])-17343.0*phiLx[0]+1170960.0*phiC[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1122732.653974222*rdx2SqVolSq[0]*phiLy[1]-1113691.348758712*rdx2SqVolSq[0]*phiLx[1]-5155056.0*rdx2SqVolSq[0]*bcVals[1]+((-2972058.0*phiLy[0])-1286154.0*phiLx[0]+6835740.0*phiC[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(639329.3979374008*rdx2SqVolCu[0]*phiLy[1]-3757176.73613406*rdx2SqVolCu[0]*phiLx[1]-1.3513464e+7*rdx2SqVolCu[0]*bcVals[1]+((-1559817.0*phiLy[0])-4278411.0*phiLx[0]+1.259496e+7*phiC[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0]*phiLx[1]-5588352.0*rdx2SqVolR4[0]*bcVals[1]+(4667544.0*phiC[0]-1873368.0*phiLx[0])*rdx2SqVolR4[0])*omega-40344.0*phiC[0]*rdx2SqVolR4[1]-1170960.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*phiC[0]*rdx2SqVolR4[0]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[1] = (((17043.37994647775*rdx2SqVolCu[1]+113483.9689119128*rdx2SqVol[0]*rdx2SqVolSq[1]+161184.6481523597*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-85680.0*rdx2SqVol[0]*rdx2SqVolSq[1])-260100.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+(92496.0*rdx2SqVolCu[1]+873696.0*rdx2SqVol[0]*rdx2SqVolSq[1]+2060172.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[1]-407045.7961851466*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1411586.767152483*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0]*rho[0])*omega*volFac+((39767.88654178142*rdx2SqVolR4[1]+404338.6007729165*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+446843.1275906566*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-20239.01368644233*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-219563.4206214686*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-522469.6620015367*rdx2SqVol[0]*rdx2SqVolCu[1])-1761980.645523667*rdx2SqVolSq[0]*rdx2SqVolSq[1]-623288.8754085057*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-168555.0*rdx2SqVol[0]*rdx2SqVolCu[1])-606900.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-337365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-20295.0*rdx2SqVol[0]*rdx2SqVolCu[1])-151200.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-240705.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(40344.0*phiLy[1]-40344.0*phiC[1])*rdx2SqVolR4[1]+(413649.0*rdx2SqVol[0]*phiLy[1]-109839.0*rdx2SqVol[0]*phiLx[1]-1170960.0*rdx2SqVol[0]*phiC[1]+560727.200239118*rdx2SqVol[0]*bcVals[1]+((-170220.7572154465*phiLy[0])-110142.8429041125*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1048338.0*rdx2SqVolSq[0]*phiLy[1]-1081578.0*rdx2SqVolSq[0]*phiLx[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]+3464794.435460782*rdx2SqVolSq[0]*bcVals[1]+((-615692.1005665087*phiLy[0])-1116705.117163882*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(464373.0*rdx2SqVolCu[0]*phiLy[1]-2599263.0*rdx2SqVolCu[0]*phiLx[1]-1.259496e+7*rdx2SqVolCu[0]*phiC[1]+6136988.564971584*rdx2SqVolCu[0]*bcVals[1]+((-350599.9924172843*phiLy[0])-2717894.290068507*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0]*phiLx[1]-4667544.0*rdx2SqVolR4[0]*phiC[1]+2199843.089677078*rdx2SqVolR4[0]*bcVals[1]-1099921.544838539*phiLx[0]*rdx2SqVolR4[0])*omega+40344.0*phiC[1]*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]*phiC[1])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((56700.41523657476*rdx2SqVol[0]*rdx2SqVolSq[1]+492907.0188179509*rdx2SqVolSq[0]*rdx2SqVol[1]+1136585.596333157*rdx2SqVolCu[0])*rho[3]+((-17712.0*rdx2SqVolCu[1])-472896.0*rdx2SqVol[0]*rdx2SqVolSq[1]-2197476.0*rdx2SqVolSq[0]*rdx2SqVol[1]-2773008.0*rdx2SqVolCu[0])*rho[2]+(197904.0*rdx2SqVol[0]*rdx2SqVolSq[1]+600780.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-21588.28126553848*rho[0]*rdx2SqVolCu[1]-482840.3395243608*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-1233559.656947324*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((72862.18132199996*rdx2SqVol[0]*rdx2SqVolCu[1]+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]-686687.1311179493*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-1917.380243978746*rdx2SqVol[0]*rdx2SqVolCu[1])-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-823210.8398721432*rdx2SqVolCu[0]*rdx2SqVol[1]-1649882.317257809*rdx2SqVolR4[0])*phiLx[3]+((-45449.01319060734*rdx2SqVolR4[1])-1119902.482954654*rdx2SqVol[0]*rdx2SqVolCu[1]-4193890.830602055*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3735659.468951633*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-984.0*rdx2SqVolR4[1])+24363.0*rdx2SqVol[0]*rdx2SqVolCu[1]+659958.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1675359.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-3321.0*rdx2SqVol[0]*rdx2SqVolCu[1])-156186.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-998973.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1873368.0*rdx2SqVolR4[0])*phiLx[2]+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[2]+(72447.0*rdx2SqVol[0]*phiLy[1]-2337.0*rdx2SqVol[0]*phiLx[1]-97147.26569492316*rdx2SqVol[0]*bcVals[1]+(52621.43558475006*phiLy[0]-4047.802737288466*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-287280.0*rdx2SqVolSq[0]*phiLx[1])-973052.2872857346*rdx2SqVolSq[0]*bcVals[1]+(812719.8081306986*phiLy[0]-326193.6644878315*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-779247.0*rdx2SqVolCu[0]*phiLy[1])-846963.0*rdx2SqVolCu[0]*phiLx[1]-1908983.261663653*rdx2SqVolCu[0]*bcVals[1]+(1901183.83687717*phiLy[0]-946692.2060453439*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[2]))/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 
  phiC[3] = (((17712.0*rdx2SqVolCu[1]+375744.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1352916.0*rdx2SqVolSq[0]*rdx2SqVol[1]+825552.0*rdx2SqVolCu[0])*rho[3]+((-31093.77609747648*rdx2SqVol[0]*rdx2SqVolSq[1])-270303.8490291989*rdx2SqVolSq[0]*rdx2SqVol[1]-623288.8754085057*rdx2SqVolCu[0])*rho[2]+(21588.28126553848*rdx2SqVolCu[1]+143746.3606217562*rdx2SqVol[0]*rdx2SqVolSq[1]+204167.220992989*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-108528.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-329460.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((984.0*rdx2SqVolR4[1]-149207.0*rdx2SqVol[0]*rdx2SqVolCu[1]-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-498771.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-21033.0*rdx2SqVol[0]*rdx2SqVolCu[1])-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1635849.0*rdx2SqVolCu[0]*rdx2SqVol[1]-1047816.0*rdx2SqVolR4[0])*phiLx[3]+((-40344.0*rdx2SqVolR4[1])-1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-4667544.0*rdx2SqVolR4[0])*phiC[3]+((-150416.0*rdx2SqVol[0]*rdx2SqVolCu[1])-693600.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-839664.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-39956.68007980643*rdx2SqVol[0]*rdx2SqVolCu[1])-15016.88050162216*rdx2SqVolSq[0]*rdx2SqVolSq[1]+376570.3622259722*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-21091.18268376621*rdx2SqVol[0]*rdx2SqVolCu[1])-453260.3758326994*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1661713.956324312*rdx2SqVolCu[0]*rdx2SqVol[1]-1099921.544838539*rdx2SqVolR4[0])*phiLx[2]+((-176754.0528615963*rdx2SqVol[0]*phiLy[1])-25636.08400282695*rdx2SqVol[0]*phiLx[1]+130872.0*rdx2SqVol[0]*bcVals[1]+((-39729.0*phiLy[0])-25707.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-812719.8081306986*rdx2SqVolSq[0]*phiLy[1])-182447.3038660752*rdx2SqVolSq[0]*phiLx[1]+383040.0*rdx2SqVolSq[0]*bcVals[1]-191520.0*phiLx[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-566001.2949481651*rdx2SqVolCu[0]*phiLy[1])-278113.666120527*rdx2SqVolCu[0]*phiLx[1]-244872.0*rdx2SqVolCu[0]*bcVals[1]+(427329.0*phiLy[0]-304893.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0])*phiC[3])/(40344.0*rdx2SqVolR4[1]+1170960.0*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.259496e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+4667544.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((474300.0*rdx2SqVol[0]*rdx2SqVolSq[1]+156240.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1136585.596333157*rdx2SqVolCu[1]+2574069.987160411*rdx2SqVol[0]*rdx2SqVolSq[1]+742259.9812787967*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-973862.8870636768*rdx2SqVol[0]*rdx2SqVolSq[1])-381189.7417297584*rdx2SqVolSq[0]*rdx2SqVol[1]-17043.37994647775*rdx2SqVolCu[0])*rho[1]-2773008.0*rho[0]*rdx2SqVolCu[1]-6470652.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-2145504.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-668655.0*rdx2SqVol[0]*rdx2SqVolCu[1])-226800.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1845.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(615195.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1106700.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+307365.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-5588352.0*rdx2SqVolR4[1])-1.3513464e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-5155056.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-416232.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1649882.317257809*rdx2SqVolR4[1])-3757176.73613406*rdx2SqVol[0]*rdx2SqVolCu[1]-1113691.348758712*rdx2SqVolSq[0]*rdx2SqVolSq[1]-10012.98571855568*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(639329.3979374008*rdx2SqVol[0]*rdx2SqVolCu[1]+1122732.653974222*rdx2SqVolSq[0]*rdx2SqVolSq[1]+310402.557275226*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(4667544.0*phiC[0]-1873368.0*phiLy[0])*rdx2SqVolR4[1]+((-747388.5837200081*rdx2SqVol[0]*phiLy[1])-1500934.608060923*rdx2SqVol[0]*phiLx[1]-2773008.0*rdx2SqVol[0]*bcVals[1]+((-4278411.0*phiLy[0])-1559817.0*phiLx[0]+1.259496e+7*phiC[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-257521.3140693406*rdx2SqVolSq[0]*phiLy[1])-2913967.637637727*rdx2SqVolSq[0]*phiLx[1]-7782592.0*rdx2SqVolSq[0]*bcVals[1]+((-1286154.0*phiLy[0])-2972058.0*phiLx[0]+6835740.0*phiC[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-3195.633739964578*rdx2SqVolCu[0]*phiLy[1])-930985.9693223094*rdx2SqVolCu[0]*phiLx[1]-2659024.0*rdx2SqVolCu[0]*bcVals[1]+((-17343.0*phiLy[0])-945501.0*phiLx[0]+1170960.0*phiC[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolR4[0]*phiLx[1]-115456.0*rdx2SqVolR4[0]*bcVals[1]+(40344.0*phiC[0]-40344.0*phiLx[0])*rdx2SqVolR4[0])*omega-4667544.0*phiC[0]*rdx2SqVolR4[1]-1.259496e+7*phiC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*phiC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*phiC[0]*rdx2SqVolR4[0]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((1136585.596333157*rdx2SqVolCu[1]+492907.0188179509*rdx2SqVol[0]*rdx2SqVolSq[1]+56700.41523657476*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(600780.0*rdx2SqVol[0]*rdx2SqVolSq[1]+197904.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-2773008.0*rdx2SqVolCu[1])-2197476.0*rdx2SqVol[0]*rdx2SqVolSq[1]-472896.0*rdx2SqVolSq[0]*rdx2SqVol[1]-17712.0*rdx2SqVolCu[0])*rho[1]-1233559.656947324*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-482840.3395243608*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-21588.28126553848*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-1649882.317257809*rdx2SqVolR4[1])-823210.8398721432*rdx2SqVol[0]*rdx2SqVolCu[1]-118524.2367619383*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1917.380243978746*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-686687.1311179493*rdx2SqVol[0]*rdx2SqVolCu[1])+27383.72326766395*rdx2SqVolSq[0]*rdx2SqVolSq[1]+72862.18132199996*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-1908983.261663653*rdx2SqVol[0]*rdx2SqVolCu[1])-973052.2872857346*rdx2SqVolSq[0]*rdx2SqVolSq[1]-97147.26569492316*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-846963.0*rdx2SqVol[0]*rdx2SqVolCu[1])-287280.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2337.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(72447.0*rdx2SqVolCu[0]*rdx2SqVol[1]-779247.0*rdx2SqVol[0]*rdx2SqVolCu[1])*phiLx[2]+(4667544.0*phiC[1]-1873368.0*phiLy[1])*rdx2SqVolR4[1]+((-998973.0*rdx2SqVol[0]*phiLy[1])+1675359.0*rdx2SqVol[0]*phiLx[1]+1.259496e+7*rdx2SqVol[0]*phiC[1]-3735659.468951633*rdx2SqVol[0]*bcVals[1]+(1901183.83687717*phiLx[0]-946692.2060453439*phiLy[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-156186.0*rdx2SqVolSq[0]*phiLy[1])+659958.0*rdx2SqVolSq[0]*phiLx[1]+6835740.0*rdx2SqVolSq[0]*phiC[1]-4193890.830602055*rdx2SqVolSq[0]*bcVals[1]+(812719.8081306986*phiLx[0]-326193.6644878315*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-3321.0*rdx2SqVolCu[0]*phiLy[1])+24363.0*rdx2SqVolCu[0]*phiLx[1]+1170960.0*rdx2SqVolCu[0]*phiC[1]-1119902.482954654*rdx2SqVolCu[0]*bcVals[1]+(52621.43558475006*phiLx[0]-4047.802737288466*phiLy[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-984.0*rdx2SqVolR4[0]*phiLx[1]+40344.0*rdx2SqVolR4[0]*phiC[1]-45449.01319060734*rdx2SqVolR4[0]*bcVals[1])*omega-4667544.0*phiC[1]*rdx2SqVolR4[1]-1.259496e+7*rdx2SqVol[0]*phiC[1]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*phiC[1]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0]*phiC[1]))/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[2] = (((161184.6481523597*rdx2SqVol[0]*rdx2SqVolSq[1]+113483.9689119128*rdx2SqVolSq[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolCu[0])*rho[3]+(825552.0*rdx2SqVolCu[1]+2060172.0*rdx2SqVol[0]*rdx2SqVolSq[1]+873696.0*rdx2SqVolSq[0]*rdx2SqVol[1]+92496.0*rdx2SqVolCu[0])*rho[2]+((-260100.0*rdx2SqVol[0]*rdx2SqVolSq[1])-85680.0*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-623288.8754085057*rho[0]*rdx2SqVolCu[1]-1411586.767152483*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-407045.7961851466*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-219563.4206214686*rdx2SqVol[0]*rdx2SqVolCu[1])-144037.3451574278*rdx2SqVolSq[0]*rdx2SqVolSq[1]-20239.01368644233*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(446843.1275906566*rdx2SqVol[0]*rdx2SqVolCu[1]+1017718.413511321*rdx2SqVolSq[0]*rdx2SqVolSq[1]+404338.6007729165*rdx2SqVolCu[0]*rdx2SqVol[1]+39767.88654178142*rdx2SqVolR4[0])*phiLx[3]+(2199843.089677078*rdx2SqVolR4[1]+6136988.564971584*rdx2SqVol[0]*rdx2SqVolCu[1]+3464794.435460782*rdx2SqVolSq[0]*rdx2SqVolSq[1]+560727.200239118*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1047816.0*rdx2SqVolR4[1])-2599263.0*rdx2SqVol[0]*rdx2SqVolCu[1]-1081578.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-109839.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(464373.0*rdx2SqVol[0]*rdx2SqVolCu[1]+1048338.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+413649.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiLx[2]+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[2]-1099921.544838539*phiLy[0]*rdx2SqVolR4[1]+((-240705.0*rdx2SqVol[0]*phiLy[1])-337365.0*rdx2SqVol[0]*phiLx[1]-623288.8754085057*rdx2SqVol[0]*bcVals[1]+((-2717894.290068507*phiLy[0])-350599.9924172843*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-151200.0*rdx2SqVolSq[0]*phiLy[1])-606900.0*rdx2SqVolSq[0]*phiLx[1]-1761980.645523667*rdx2SqVolSq[0]*bcVals[1]+((-1116705.117163882*phiLy[0])-615692.1005665087*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-20295.0*rdx2SqVolCu[0]*phiLy[1])-168555.0*rdx2SqVolCu[0]*phiLx[1]-522469.6620015367*rdx2SqVolCu[0]*bcVals[1]+((-110142.8429041125*phiLy[0])-170220.7572154465*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[2])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 
  phiC[3] = (((825552.0*rdx2SqVolCu[1]+1352916.0*rdx2SqVol[0]*rdx2SqVolSq[1]+375744.0*rdx2SqVolSq[0]*rdx2SqVol[1]+17712.0*rdx2SqVolCu[0])*rho[3]+(204167.220992989*rdx2SqVol[0]*rdx2SqVolSq[1]+143746.3606217562*rdx2SqVolSq[0]*rdx2SqVol[1]+21588.28126553848*rdx2SqVolCu[0])*rho[2]+((-623288.8754085057*rdx2SqVolCu[1])-270303.8490291989*rdx2SqVol[0]*rdx2SqVolSq[1]-31093.77609747648*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-329460.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-108528.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-1047816.0*rdx2SqVolR4[1])-1635849.0*rdx2SqVol[0]*rdx2SqVolCu[1]-449562.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-21033.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-498771.0*rdx2SqVol[0]*rdx2SqVolCu[1])-706878.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-149207.0*rdx2SqVolCu[0]*rdx2SqVol[1]+984.0*rdx2SqVolR4[0])*phiLx[3]+((-4667544.0*rdx2SqVolR4[1])-1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]-40344.0*rdx2SqVolR4[0])*phiC[3]+((-244872.0*rdx2SqVol[0]*rdx2SqVolCu[1])+383040.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+130872.0*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-278113.666120527*rdx2SqVol[0]*rdx2SqVolCu[1])-182447.3038660752*rdx2SqVolSq[0]*rdx2SqVolSq[1]-25636.08400282695*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+((-566001.2949481651*rdx2SqVol[0]*rdx2SqVolCu[1])-812719.8081306986*rdx2SqVolSq[0]*rdx2SqVolSq[1]-176754.0528615963*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]-1099921.544838539*phiLy[1]*rdx2SqVolR4[1]+((-1661713.956324312*rdx2SqVol[0]*phiLy[1])+376570.3622259722*rdx2SqVol[0]*phiLx[1]-839664.0*rdx2SqVol[0]*bcVals[1]+(427329.0*phiLx[0]-304893.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-453260.3758326994*rdx2SqVolSq[0]*phiLy[1])-15016.88050162216*rdx2SqVolSq[0]*phiLx[1]-693600.0*rdx2SqVolSq[0]*bcVals[1]-191520.0*phiLy[0]*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-21091.18268376621*rdx2SqVolCu[0]*phiLy[1])-39956.68007980643*rdx2SqVolCu[0]*phiLx[1]-150416.0*rdx2SqVolCu[0]*bcVals[1]+((-25707.0*phiLy[0])-39729.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0])*phiC[3])/(4667544.0*rdx2SqVolR4[1]+1.259496e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+6835740.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1170960.0*rdx2SqVolCu[0]*rdx2SqVol[1]+40344.0*rdx2SqVolR4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(40344.0*phiLy[0]-40344.0*phiC[0])*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiLy[1]+50064.92859277839*rdx2SqVol[0]*phiLx[1]+92496.0*rdx2SqVol[0]*bcVals[1]+(293355.0*phiLy[0]+52029.0*phiLx[0]-345384.0*phiC[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+288932.0554646022*rdx2SqVolSq[0]*phiLx[1]+828720.0*rdx2SqVolSq[0]*bcVals[1]+(52029.0*phiLy[0]+293355.0*phiLx[0]-345384.0*phiC[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiLx[1]+115456.0*rdx2SqVolCu[0]*bcVals[1]+(40344.0*phiLx[0]-40344.0*phiC[0])*rdx2SqVolCu[0])*omega+40344.0*phiC[0]*rdx2SqVolCu[1]+345384.0*phiC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*phiC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*phiC[0]*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+(121032.0*phiLy[1]-121032.0*phiC[1])*rdx2SqVolCu[1]+(221031.0*rdx2SqVol[0]*phiLy[1]-167649.0*rdx2SqVol[0]*phiLx[1]-1036152.0*rdx2SqVol[0]*phiC[1]+373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLy[0]-190246.7286525579*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(29889.0*rdx2SqVolSq[0]*phiLy[1]-11367.0*rdx2SqVolSq[0]*phiLx[1]-1036152.0*rdx2SqVolSq[0]*phiC[1]+1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLy[0]-36430.22463559618*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiLx[1]-121032.0*rdx2SqVolCu[0]*phiC[1]+136347.039571822*rdx2SqVolCu[0]*bcVals[1])*omega+121032.0*phiC[1]*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*phiC[1]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*phiC[1]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]*phiC[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[2]+(35055.0*rdx2SqVol[0]*phiLy[1]+35055.0*rdx2SqVol[0]*phiLx[1]+64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-35055.0*rdx2SqVolSq[0]*phiLy[1])+188385.0*rdx2SqVolSq[0]*phiLx[1]+583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[2])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiC[3]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-190246.7286525579*rdx2SqVol[0]*phiLy[1])-39128.75979378851*rdx2SqVol[0]*phiLx[1]+87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLy[0]-44403.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-36430.22463559618*rdx2SqVolSq[0]*phiLy[1])+44657.46597154836*rdx2SqVolSq[0]*phiLx[1]+168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0])*rdx2SqVol[0])/(18.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[1] = (16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])/(18.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 
  phiC[2] = (16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]+8.660254037844386*rdx2SqVol[0]*phiLx[3]-7.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]-7.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1])/(46.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[3] = (16.0*rho[3]*volFac-7.0*rdx2SqVol[1]*phiUy[3]-7.0*rdx2SqVol[0]*phiUx[3]-7.0*rdx2SqVol[1]*phiLy[3]-7.0*rdx2SqVol[0]*phiLx[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1])/(46.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 

}

void MGpoissonJacobi2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[0]*rho[1]+288.0*rho[0]*rdx2SqVol[1]+2096.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(162.0*phiUy[0]+162.0*phiLy[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]-31.17691453623978*rdx2SqVol[0]*phiUx[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]+(1179.0*phiUy[0]+54.0*phiUx[0]+1179.0*phiLy[0]+1296.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiUx[1]+(1416.0*phiUx[0]+4224.0*bcVals[0])*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]+471.1178196587346*rdx2SqVol[0]*rho[0])*volFac+((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]+162.0*phiLy[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]-342.0*rdx2SqVol[0]*phiUx[1]+351.0*rdx2SqVol[0]*phiLy[1]+(265.0037735580381*phiUy[0]+342.9460598986376*phiUx[0]+265.0037735580381*phiLy[0]-1745.907214029428*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiUx[1]+(831.384387633061*phiUx[0]-1662.768775266122*bcVals[0])*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = ((859.0972005541631*rdx2SqVol[0]*rho[3]+(736.0*rdx2SqVol[1]+2096.0*rdx2SqVol[0])*rho[2])*volFac-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(138.0*rdx2SqVol[0]*rdx2SqVol[1]+1416.0*rdx2SqVolSq[0])*phiUx[2]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(398.3716857408418*phiUy[0]-398.3716857408418*phiLy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiUy[0]-1134.493278957615*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[0]*rho[2])*volFac+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiUx[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+(876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1]+831.384387633061*rdx2SqVolSq[0])*phiUx[2]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiUy[0]-255.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[1]-864.0*rho[0]*rdx2SqVol[1]-2256.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-486.0*phiUy[0])-486.0*phiLy[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+467.6537180435967*rdx2SqVol[0]*phiUx[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+((-1269.0*phiUy[0])-486.0*phiUx[0]-1269.0*phiLy[0]+864.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiUx[1]+(2816.0*bcVals[0]-984.0*phiUx[0])*rdx2SqVolSq[0]))/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]-526.5434455009387*rdx2SqVol[0]*rho[0])*volFac+((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]+486.0*phiLy[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]-522.0*rdx2SqVol[0]*phiUx[1]+243.0*rdx2SqVol[0]*phiLy[1]+((-296.1806880942779*phiUy[0])+592.3613761885558*phiUx[0]-296.1806880942779*phiLy[0]+1163.938142686285*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiUx[1]+1108.512516844081*bcVals[0]*rdx2SqVolSq[0])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[3]+((-2208.0*rdx2SqVol[1])-2256.0*rdx2SqVol[0])*rho[2])*volFac-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-1242.0*rdx2SqVol[0]*rdx2SqVol[1])-984.0*rdx2SqVolSq[0])*phiUx[2]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1195.115057222525*phiLy[0]-1195.115057222525*phiUy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiLy[0]-1221.095819336058*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1]))/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[0]*rho[2])*volFac+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUx[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiLy[0]-285.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[1]-288.0*rho[0]*rdx2SqVol[1]-2096.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-162.0*phiUy[0])-162.0*phiLy[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]-31.17691453623978*rdx2SqVol[0]*phiLx[1]-1296.0*rdx2SqVol[0]*bcVals[1]+((-1179.0*phiUy[0])-1179.0*phiLy[0]-54.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiLx[1]-4224.0*rdx2SqVolSq[0]*bcVals[1]-1416.0*phiLx[0]*rdx2SqVolSq[0]))/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]-471.1178196587346*rdx2SqVol[0]*rho[0])*volFac+((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]+162.0*phiLy[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]+351.0*rdx2SqVol[0]*phiLy[1]-342.0*rdx2SqVol[0]*phiLx[1]+1745.907214029428*rdx2SqVol[0]*bcVals[1]+((-265.0037735580381*phiUy[0])-265.0037735580381*phiLy[0]-342.9460598986376*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiLx[1]+1662.768775266122*rdx2SqVolSq[0]*bcVals[1]-831.384387633061*phiLx[0]*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[3]+((-736.0*rdx2SqVol[1])-2096.0*rdx2SqVol[0])*rho[2])*volFac-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiLx[3]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-138.0*rdx2SqVol[0]*rdx2SqVol[1])-1416.0*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiLy[0]-398.3716857408418*phiUy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiLy[0]-1134.493278957615*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1]))/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[0]*rho[2])*volFac+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiLx[3]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+((-876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1])-831.384387633061*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiLy[0]-255.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[0]*rho[1]+864.0*rho[0]*rdx2SqVol[1]+2256.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(486.0*phiUy[0]+486.0*phiLy[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+467.6537180435967*rdx2SqVol[0]*phiLx[1]+864.0*rdx2SqVol[0]*bcVals[1]+(1269.0*phiUy[0]+1269.0*phiLy[0]+486.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiLx[1]+2816.0*rdx2SqVolSq[0]*bcVals[1]+984.0*phiLx[0]*rdx2SqVolSq[0])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]+526.5434455009387*rdx2SqVol[0]*rho[0])*volFac+((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]+486.0*phiLy[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]+243.0*rdx2SqVol[0]*phiLy[1]-522.0*rdx2SqVol[0]*phiLx[1]+1163.938142686285*rdx2SqVol[0]*bcVals[1]+(296.1806880942779*phiUy[0]+296.1806880942779*phiLy[0]-592.3613761885558*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiLx[1]+1108.512516844081*rdx2SqVolSq[0]*bcVals[1])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = ((415.6921938165305*rdx2SqVol[0]*rho[3]+(2208.0*rdx2SqVol[1]+2256.0*rdx2SqVol[0])*rho[2])*volFac-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiLx[3]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1242.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiLx[2]+(1195.115057222525*phiUy[0]-1195.115057222525*phiLy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiUy[0]-1221.095819336058*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[0]*rho[2])*volFac+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLx[3]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]-1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiUy[0]-285.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[1]*rho[2]+2096.0*rho[0]*rdx2SqVol[1]+288.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(4224.0*rdx2SqVolSq[1]+1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+1416.0*phiUy[0]*rdx2SqVolSq[1]+((-1134.493278957615*rdx2SqVol[0]*phiUx[1])+1134.493278957615*rdx2SqVol[0]*phiLx[1]+(54.0*phiUy[0]+1179.0*phiUx[0]+1179.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-155.8845726811989*rdx2SqVolSq[0]*phiUx[1]+155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+(162.0*phiUx[0]+162.0*phiLx[0])*rdx2SqVolSq[0])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = ((859.0972005541631*rdx2SqVol[1]*rho[3]+(2096.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[1])*volFac+((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+1416.0*phiUy[1]*rdx2SqVolSq[1]+(138.0*rdx2SqVol[0]*phiUy[1]-917.0*rdx2SqVol[0]*phiUx[1]-917.0*rdx2SqVol[0]*phiLx[1]+(1134.493278957615*phiUx[0]-1134.493278957615*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-322.0*rdx2SqVolSq[0]*phiUx[1]-322.0*rdx2SqVolSq[0]*phiLx[1]+(398.3716857408418*phiUx[0]-398.3716857408418*phiLx[0])*rdx2SqVolSq[0])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]+471.1178196587346*rho[0]*rdx2SqVol[1])*volFac+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]+((-1662.768775266122*rdx2SqVolSq[1])-1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+831.384387633061*phiUy[0]*rdx2SqVolSq[1]+((-255.0*rdx2SqVol[0]*phiUx[1])+255.0*rdx2SqVol[0]*phiLx[1]+(342.9460598986376*phiUy[0]+265.0037735580381*phiUx[0]+265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[1]*rho[1])*volFac+((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]+831.384387633061*phiUy[1]*rdx2SqVolSq[1]+(876.4177086298519*rdx2SqVol[0]*phiUy[1]-206.1140461006964*rdx2SqVol[0]*phiUx[1]-206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiUx[0]-255.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[2]-2256.0*rho[0]*rdx2SqVol[1]-864.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]-984.0*phiUy[0]*rdx2SqVolSq[1]+(1221.095819336058*rdx2SqVol[0]*phiUx[1]-1221.095819336058*rdx2SqVol[0]*phiLx[1]+((-486.0*phiUy[0])-1269.0*phiUx[0]-1269.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0]*phiUx[1]-467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+((-486.0*phiUx[0])-486.0*phiLx[0])*rdx2SqVolSq[0]))/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[3]+((-2256.0*rdx2SqVol[1])-2208.0*rdx2SqVol[0])*rho[1])*volFac+(969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]-984.0*phiUy[1]*rdx2SqVolSq[1]+((-1242.0*rdx2SqVol[0]*phiUy[1])+987.0*rdx2SqVol[0]*phiUx[1]+987.0*rdx2SqVol[0]*phiLx[1]+(1221.095819336058*phiLx[0]-1221.095819336058*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+966.0*rdx2SqVolSq[0]*phiUx[1]+966.0*rdx2SqVolSq[0]*phiLx[1]+(1195.115057222525*phiLx[0]-1195.115057222525*phiUx[0])*rdx2SqVolSq[0]))/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]-526.5434455009387*rho[0]*rdx2SqVol[1])*volFac+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(285.0*rdx2SqVol[0]*phiUx[1]-285.0*rdx2SqVol[0]*phiLx[1]+(592.3613761885558*phiUy[0]-296.1806880942779*phiUx[0]-296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[1]*rho[1])*volFac+(24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+(1513.812405815199*rdx2SqVol[0]*phiUy[1]+230.3627574066607*rdx2SqVol[0]*phiUx[1]+230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiLx[0]-285.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[2]-2096.0*rho[0]*rdx2SqVol[1]-288.0*rdx2SqVol[0]*rho[0])*volFac-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-4224.0*rdx2SqVolSq[1])-1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]-1416.0*phiLy[0]*rdx2SqVolSq[1]+(1134.493278957615*rdx2SqVol[0]*phiUx[1]-1134.493278957615*rdx2SqVol[0]*phiLx[1]+((-1179.0*phiUx[0])-54.0*phiLy[0]-1179.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0]*phiUx[1]-155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+((-162.0*phiUx[0])-162.0*phiLx[0])*rdx2SqVolSq[0]))/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[3]+((-2096.0*rdx2SqVol[1])-736.0*rdx2SqVol[0])*rho[1])*volFac-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]-1416.0*phiLy[1]*rdx2SqVolSq[1]+(917.0*rdx2SqVol[0]*phiUx[1]-138.0*rdx2SqVol[0]*phiLy[1]+917.0*rdx2SqVol[0]*phiLx[1]+(1134.493278957615*phiLx[0]-1134.493278957615*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+322.0*rdx2SqVolSq[0]*phiUx[1]+322.0*rdx2SqVolSq[0]*phiLx[1]+(398.3716857408418*phiLx[0]-398.3716857408418*phiUx[0])*rdx2SqVolSq[0]))/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]-471.1178196587346*rho[0]*rdx2SqVol[1])*volFac+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+(1662.768775266122*rdx2SqVolSq[1]+1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[0]*rdx2SqVolSq[1]+(255.0*rdx2SqVol[0]*phiUx[1]-255.0*rdx2SqVol[0]*phiLx[1]+((-265.0037735580381*phiUx[0])-342.9460598986376*phiLy[0]-265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[1]*rho[1])*volFac+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[1]*rdx2SqVolSq[1]+(206.1140461006964*rdx2SqVol[0]*phiUx[1]-876.4177086298519*rdx2SqVol[0]*phiLy[1]+206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiLx[0]-255.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[1]*rho[2]+2256.0*rho[0]*rdx2SqVol[1]+864.0*rdx2SqVol[0]*rho[0])*volFac-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+984.0*phiLy[0]*rdx2SqVolSq[1]+((-1221.095819336058*rdx2SqVol[0]*phiUx[1])+1221.095819336058*rdx2SqVol[0]*phiLx[1]+(1269.0*phiUx[0]+486.0*phiLy[0]+1269.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-467.6537180435967*rdx2SqVolSq[0]*phiUx[1]+467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+(486.0*phiUx[0]+486.0*phiLx[0])*rdx2SqVolSq[0])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = ((415.6921938165305*rdx2SqVol[1]*rho[3]+(2256.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[1])*volFac-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]+(969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+984.0*phiLy[1]*rdx2SqVolSq[1]+((-987.0*rdx2SqVol[0]*phiUx[1])+1242.0*rdx2SqVol[0]*phiLy[1]-987.0*rdx2SqVol[0]*phiLx[1]+(1221.095819336058*phiUx[0]-1221.095819336058*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-966.0*rdx2SqVolSq[0]*phiUx[1]-966.0*rdx2SqVolSq[0]*phiLx[1]+(1195.115057222525*phiUx[0]-1195.115057222525*phiLx[0])*rdx2SqVolSq[0])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]+526.5434455009387*rho[0]*rdx2SqVol[1])*volFac+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+((-285.0*rdx2SqVol[0]*phiUx[1])+285.0*rdx2SqVol[0]*phiLx[1]+(296.1806880942779*phiUx[0]-592.3613761885558*phiLy[0]+296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[1]*rho[1])*volFac+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+(24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+((-230.3627574066607*rdx2SqVol[0]*phiUx[1])-1513.812405815199*rdx2SqVol[0]*phiLy[1]-230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiUx[0]-285.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonJacobi2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+624456.0*phiUy[0]*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiUx[1]+(3368931.0*phiUy[0]+173313.0*phiUx[0]+4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiUx[1]+(173313.0*phiUy[0]+3368931.0*phiUx[0]+1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+(624456.0*phiUx[0]+1862784.0*bcVals[0])*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+624456.0*phiUy[1]*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiUx[1]+(847040.3948826761*phiUy[0]+1100685.379244677*phiUx[0]-5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiUx[1]+(38955.5547130316*phiUy[0]+2275410.734360502*phiUx[0]-5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiUx[1]+(366640.5149461798*phiUx[0]-733281.0298923596*bcVals[0])*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiUx[1]+(2275410.734360502*phiUy[0]+38955.5547130316*phiUx[0]+934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiUx[1]+(1100685.379244677*phiUy[0]+847040.3948826761*phiUx[0]+1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiUx[1]+(90933.0*phiUy[0]+82467.0*phiUx[0]-419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiUx[1]+(82467.0*phiUy[0]+90933.0*phiUx[0]+73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1.1265816e+7*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(1.064828809836548e+7*rdx2SqVolSq[1]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(1.064828809836548e+7*rdx2SqVolSq[0]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-9250416.0*rho[0]*rdx2SqVolSq[1]+1.7580136e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-9250416.0*rdx2SqVolSq[0]*rho[0])*volFac+(8660259.0*rdx2SqVol[0]*rdx2SqVolSq[1]-7080939.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(8660259.0*rdx2SqVolSq[0]*rdx2SqVol[1]-7080939.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(1492750.66799516*rdx2SqVolCu[1]-2750120.827394134*rdx2SqVol[0]*rdx2SqVolSq[1]+6151376.711030058*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(7130550.065869739*rdx2SqVol[0]*rdx2SqVolSq[1]-7586460.479438022*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.0336768e+7*rdx2SqVolCu[1]-5.7997776e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1893392e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-430920.0*phiUy[0]*rdx2SqVolCu[1]+((-7586460.479438022*rdx2SqVol[0]*phiUy[1])+6151376.711030058*rdx2SqVol[0]*phiUx[1]+(711555.0*phiUy[0]-6194475.0*phiUx[0]+1.1893392e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(7130550.065869739*rdx2SqVolSq[0]*phiUy[1]-2750120.827394134*rdx2SqVolSq[0]*phiUx[1]+((-6194475.0*phiUy[0])+711555.0*phiUx[0]-5.7997776e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1492750.66799516*rdx2SqVolCu[0]*phiUx[1]+(3.0336768e+7*bcVals[0]-430920.0*phiUx[0])*rdx2SqVolCu[0])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[1] = (((1.064828809836548e+7*rdx2SqVolSq[1]-1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-9250416.0*rdx2SqVolSq[1])+158424.0*rdx2SqVol[0]*rdx2SqVol[1]-172368.0*rdx2SqVolSq[0])*rho[1]-1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+9056020.719170641*rdx2SqVolSq[0]*rho[0])*volFac+(1492750.66799516*rdx2SqVolCu[1]+4633060.973115181*rdx2SqVol[0]*rdx2SqVolSq[1]+114621.9262924855*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-4658626.043034897*rdx2SqVol[0]*rdx2SqVolSq[1])-1632937.664207363*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(7365267.0*rdx2SqVol[0]*rdx2SqVolSq[1]-6022107.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(5763555.0*rdx2SqVol[0]*rdx2SqVolSq[1]+553725.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.894013253264102e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.164345521036225e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-430920.0*phiUy[1]*rdx2SqVolCu[1]+((-5756175.0*rdx2SqVol[0]*phiUy[1])+4047057.0*rdx2SqVol[0]*phiUx[1]+((-6452036.482512711*phiUy[0])-5006934.532233768*phiUx[0]-1.602219050314806e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-115425.0*rdx2SqVolSq[0]*phiUy[1])+1.1487735e+7*rdx2SqVolSq[0]*phiUx[1]+(6064299.588730339*phiUy[0]-1.155226793149618e+7*phiUx[0]+2.261939189589393e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-5946696.0*rdx2SqVolCu[0]*phiUx[1]+(5971002.671980642*phiUx[0]-1.194200534396128e+7*bcVals[0])*rdx2SqVolCu[0])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1]-1.064828809836548e+7*rdx2SqVolSq[0])*rho[3]+(172368.0*rdx2SqVolSq[1]-158424.0*rdx2SqVol[0]*rdx2SqVol[1]+9250416.0*rdx2SqVolSq[0])*rho[2]-9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-9056020.719170641*rho[0]*rdx2SqVolSq[1]+1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(1632937.664207363*rdx2SqVol[0]*rdx2SqVolSq[1]+4658626.043034897*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-114621.9262924855*rdx2SqVol[0]*rdx2SqVolSq[1])-4633060.973115181*rdx2SqVolSq[0]*rdx2SqVol[1]-1492750.66799516*rdx2SqVolCu[0])*phiUx[3]+(5946696.0*rdx2SqVolCu[1]-1.1487735e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4047057.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(115425.0*rdx2SqVol[0]*rdx2SqVolSq[1]+5756175.0*rdx2SqVolSq[0]*rdx2SqVol[1]+430920.0*rdx2SqVolCu[0])*phiUx[2]+(1.194200534396128e+7*rdx2SqVolCu[1]-2.261939189589393e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.602219050314806e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-5971002.671980642*phiUy[0]*rdx2SqVolCu[1]+((-553725.0*rdx2SqVol[0]*phiUy[1])+6022107.0*rdx2SqVol[0]*phiUx[1]+(1.155226793149618e+7*phiUy[0]-6064299.588730339*phiUx[0]+1.164345521036225e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5763555.0*rdx2SqVolSq[0]*phiUy[1])-7365267.0*rdx2SqVolSq[0]*phiUx[1]+(5006934.532233768*phiUy[0]+6452036.482512711*phiUx[0]-3.894013253264102e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[3] = -(1.0*(((57456.0*rdx2SqVolSq[1]+3025032.0*rdx2SqVol[0]*rdx2SqVol[1]+57456.0*rdx2SqVolSq[0])*rho[3]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[0])*rho[2]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[1])*rho[1]-2716168.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(1982232.0*rdx2SqVolCu[1]-3365199.0*rdx2SqVol[0]*rdx2SqVolSq[1]-25137.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-25137.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3365199.0*rdx2SqVolSq[0]*rdx2SqVol[1]+1982232.0*rdx2SqVolCu[0])*phiUx[3]+(462920.0231865111*rdx2SqVol[0]*rdx2SqVolSq[1]+1320669.688212385*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(31098.97224989918*rdx2SqVol[0]*rdx2SqVolSq[1]+3693399.16129776*rdx2SqVolSq[0]*rdx2SqVol[1]-1990334.223993547*rdx2SqVolCu[0])*phiUx[2]+(8810256.0*rdx2SqVol[0]*rdx2SqVolSq[1]-5228496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-1990334.223993547*phiUy[1]*rdx2SqVolCu[1]+(3693399.16129776*rdx2SqVol[0]*phiUy[1]+1320669.688212385*rdx2SqVol[0]*phiUx[1]+((-156975.0*phiUy[0])-1633905.0*phiUx[0]-5228496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(31098.97224989918*rdx2SqVolSq[0]*phiUy[1]+462920.0231865111*rdx2SqVolSq[0]*phiUx[1]+((-1633905.0*phiUy[0])-156975.0*phiUx[0]+8810256.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(4912488.0*rdx2SqVolCu[1]-9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]-9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]+4912488.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1.1265816e+7*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(1.064828809836548e+7*rdx2SqVolSq[1]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(1.064828809836548e+7*rdx2SqVolSq[0]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-9250416.0*rho[0]*rdx2SqVolSq[1]+1.7580136e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-9250416.0*rdx2SqVolSq[0]*rho[0])*volFac+(8660259.0*rdx2SqVol[0]*rdx2SqVolSq[1]-7080939.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(8660259.0*rdx2SqVolSq[0]*rdx2SqVol[1]-7080939.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(1492750.66799516*rdx2SqVolCu[1]-2750120.827394134*rdx2SqVol[0]*rdx2SqVolSq[1]+6151376.711030058*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(7130550.065869739*rdx2SqVol[0]*rdx2SqVolSq[1]-7586460.479438022*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.0336768e+7*rdx2SqVolCu[1]-5.7997776e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1893392e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-430920.0*phiUy[0]*rdx2SqVolCu[1]+((-7586460.479438022*rdx2SqVol[0]*phiUy[1])+6151376.711030058*rdx2SqVol[0]*phiUx[1]+(711555.0*phiUy[0]-6194475.0*phiUx[0]+1.1893392e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(7130550.065869739*rdx2SqVolSq[0]*phiUy[1]-2750120.827394134*rdx2SqVolSq[0]*phiUx[1]+((-6194475.0*phiUy[0])+711555.0*phiUx[0]-5.7997776e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1492750.66799516*rdx2SqVolCu[0]*phiUx[1]+(3.0336768e+7*bcVals[0]-430920.0*phiUx[0])*rdx2SqVolCu[0])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[1] = (((1.064828809836548e+7*rdx2SqVolSq[1]-1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-9250416.0*rdx2SqVolSq[1])+158424.0*rdx2SqVol[0]*rdx2SqVol[1]-172368.0*rdx2SqVolSq[0])*rho[1]-1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+9056020.719170641*rdx2SqVolSq[0]*rho[0])*volFac+(1492750.66799516*rdx2SqVolCu[1]+4633060.973115181*rdx2SqVol[0]*rdx2SqVolSq[1]+114621.9262924855*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-4658626.043034897*rdx2SqVol[0]*rdx2SqVolSq[1])-1632937.664207363*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(7365267.0*rdx2SqVol[0]*rdx2SqVolSq[1]-6022107.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(5763555.0*rdx2SqVol[0]*rdx2SqVolSq[1]+553725.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.894013253264102e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.164345521036225e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-430920.0*phiUy[1]*rdx2SqVolCu[1]+((-5756175.0*rdx2SqVol[0]*phiUy[1])+4047057.0*rdx2SqVol[0]*phiUx[1]+((-6452036.482512711*phiUy[0])-5006934.532233768*phiUx[0]-1.602219050314806e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-115425.0*rdx2SqVolSq[0]*phiUy[1])+1.1487735e+7*rdx2SqVolSq[0]*phiUx[1]+(6064299.588730339*phiUy[0]-1.155226793149618e+7*phiUx[0]+2.261939189589393e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-5946696.0*rdx2SqVolCu[0]*phiUx[1]+(5971002.671980642*phiUx[0]-1.194200534396128e+7*bcVals[0])*rdx2SqVolCu[0])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1]-1.064828809836548e+7*rdx2SqVolSq[0])*rho[3]+(172368.0*rdx2SqVolSq[1]-158424.0*rdx2SqVol[0]*rdx2SqVol[1]+9250416.0*rdx2SqVolSq[0])*rho[2]-9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-9056020.719170641*rho[0]*rdx2SqVolSq[1]+1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(1632937.664207363*rdx2SqVol[0]*rdx2SqVolSq[1]+4658626.043034897*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-114621.9262924855*rdx2SqVol[0]*rdx2SqVolSq[1])-4633060.973115181*rdx2SqVolSq[0]*rdx2SqVol[1]-1492750.66799516*rdx2SqVolCu[0])*phiUx[3]+(5946696.0*rdx2SqVolCu[1]-1.1487735e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4047057.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(115425.0*rdx2SqVol[0]*rdx2SqVolSq[1]+5756175.0*rdx2SqVolSq[0]*rdx2SqVol[1]+430920.0*rdx2SqVolCu[0])*phiUx[2]+(1.194200534396128e+7*rdx2SqVolCu[1]-2.261939189589393e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.602219050314806e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-5971002.671980642*phiUy[0]*rdx2SqVolCu[1]+((-553725.0*rdx2SqVol[0]*phiUy[1])+6022107.0*rdx2SqVol[0]*phiUx[1]+(1.155226793149618e+7*phiUy[0]-6064299.588730339*phiUx[0]+1.164345521036225e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5763555.0*rdx2SqVolSq[0]*phiUy[1])-7365267.0*rdx2SqVolSq[0]*phiUx[1]+(5006934.532233768*phiUy[0]+6452036.482512711*phiUx[0]-3.894013253264102e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[3] = -(1.0*(((57456.0*rdx2SqVolSq[1]+3025032.0*rdx2SqVol[0]*rdx2SqVol[1]+57456.0*rdx2SqVolSq[0])*rho[3]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[0])*rho[2]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[1])*rho[1]-2716168.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(1982232.0*rdx2SqVolCu[1]-3365199.0*rdx2SqVol[0]*rdx2SqVolSq[1]-25137.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-25137.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3365199.0*rdx2SqVolSq[0]*rdx2SqVol[1]+1982232.0*rdx2SqVolCu[0])*phiUx[3]+(462920.0231865111*rdx2SqVol[0]*rdx2SqVolSq[1]+1320669.688212385*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(31098.97224989918*rdx2SqVol[0]*rdx2SqVolSq[1]+3693399.16129776*rdx2SqVolSq[0]*rdx2SqVol[1]-1990334.223993547*rdx2SqVolCu[0])*phiUx[2]+(8810256.0*rdx2SqVol[0]*rdx2SqVolSq[1]-5228496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-1990334.223993547*phiUy[1]*rdx2SqVolCu[1]+(3693399.16129776*rdx2SqVol[0]*phiUy[1]+1320669.688212385*rdx2SqVol[0]*phiUx[1]+((-156975.0*phiUy[0])-1633905.0*phiUx[0]-5228496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(31098.97224989918*rdx2SqVolSq[0]*phiUy[1]+462920.0231865111*rdx2SqVolSq[0]*phiUx[1]+((-1633905.0*phiUy[0])-156975.0*phiUx[0]+8810256.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(4912488.0*rdx2SqVolCu[1]-9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]-9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]+4912488.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+40344.0*phiUy[0]*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiUx[1]+(293355.0*phiUy[0]+52029.0*phiUx[0]-92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+(52029.0*phiUy[0]+293355.0*phiUx[0]-828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+(40344.0*phiUx[0]-115456.0*bcVals[0])*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-121032.0*phiUy[1]*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiUx[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiUx[0]-373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiUx[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiUx[0]-1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiUx[1]-136347.039571822*bcVals[0]*rdx2SqVolCu[0]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiUx[1]+((-36430.22463559618*phiUy[0])+36430.22463559618*phiUx[0]-64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiUx[1]+((-190246.7286525579*phiUy[0])+190246.7286525579*phiUx[0]-583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUy[2]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiUx[1]+(44403.0*phiUy[0]-44403.0*phiUx[0]-87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiUx[1]+((-44403.0*phiUy[0])+44403.0*phiUx[0]-168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-624456.0*phiLy[0]*rdx2SqVolCu[1]+(100062.3072040616*rdx2SqVol[0]*phiUx[1]-1544603.07302135*rdx2SqVol[0]*phiLy[1]+((-173313.0*phiUx[0])-3368931.0*phiLy[0]-4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2951378.203030407*rdx2SqVolSq[0]*phiUx[1]-71036.59977082233*rdx2SqVolSq[0]*phiLy[1]+((-3368931.0*phiUx[0])-173313.0*phiLy[0]-1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+((-624456.0*phiUx[0])-1862784.0*bcVals[0])*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-624456.0*phiLy[1]*rdx2SqVolCu[1]+(1097649.0*rdx2SqVol[0]*phiUx[1]-722367.0*rdx2SqVol[0]*phiLy[1]+((-1100685.379244677*phiUx[0])-847040.3948826761*phiLy[0]+5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2182239.0*rdx2SqVolSq[0]*phiUx[1]-51597.0*rdx2SqVolSq[0]*phiLy[1]+((-2275410.734360502*phiUx[0])-38955.5547130316*phiLy[0]+5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiUx[1]+(733281.0298923596*bcVals[0]-366640.5149461798*phiUx[0])*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+(22491.0*rdx2SqVol[0]*phiUx[1]-497457.0*rdx2SqVol[0]*phiLy[1]+((-38955.5547130316*phiUx[0])-2275410.734360502*phiLy[0]-934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(757809.0*rdx2SqVolSq[0]*phiUx[1]-451143.0*rdx2SqVolSq[0]*phiLy[1]+((-847040.3948826761*phiUx[0])-1100685.379244677*phiLy[0]-1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[3]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+(82239.50439417786*rdx2SqVol[0]*phiUx[1]-474351.5585164657*rdx2SqVol[0]*phiLy[1]+((-82467.0*phiUx[0])-90933.0*phiLy[0]+419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(82946.18112366594*rdx2SqVolSq[0]*phiUx[1]-109228.3200777161*rdx2SqVolSq[0]*phiLy[1]+((-90933.0*phiUx[0])-82467.0*phiLy[0]-73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((6.447991200000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.9688856e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.887905122726076e+7*rdx2SqVolCu[1])-5.601520208069406e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-3.571379299595986e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.214413387187174e+9*rdx2SqVol[0]*rdx2SqVolSq[1])-4.580379954982652e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.840632657003176e+8*rdx2SqVolCu[0])*rho[1]+6.7183704e+8*rho[0]*rdx2SqVolCu[1]+1.307751256e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.97417736e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.59900048e+8*rdx2SqVolCu[0]*rho[0])*volFac+(1.9204101e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.5135219e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-6.07505661e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.513594099999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(5.178539520000001e+8*rdx2SqVolR4[1]+1.008237552e+9*rdx2SqVol[0]*rdx2SqVolCu[1]+2.03977536e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.30827312e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.933865037539783e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2609403.823634816*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.325858046406458e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(3.312597052538735e+8*rdx2SqVolR4[1]+6.446659452024169e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.629314075353864e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.395957580471022e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+3.54553416e+8*phiLy[0]*rdx2SqVolR4[1]+((-4.467607425939947e+8*rdx2SqVol[0]*phiUx[1])-6.504368232599894e+8*rdx2SqVol[0]*phiLy[1]+(4.49890875e+8*phiUx[0]+6.899945009999999e+8*phiLy[0]-8.6379048e+8*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.383275325638488e+8*rdx2SqVolSq[0]*phiUx[1])-2.41723701273902e+8*rdx2SqVolSq[0]*phiLy[1]+(2.1840576e+8*phiUx[0]+1.74784848e+8*phiLy[0]-3.78482016e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2613649.080164168*rdx2SqVolCu[0]*phiUx[1]+9.098581884049787e+7*rdx2SqVolCu[0]*phiLy[1]+((-4.8756675e+7*phiUx[0])-7.904150099999999e+7*phiLy[0]-1.175739312e+9*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+2.580326154677349e+7*rdx2SqVolR4[0]*phiUx[1]+(5.243927039999999e+8*bcVals[0]-7448760.0*phiUx[0])*rdx2SqVolR4[0])/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((2.887905122726076e+7*rdx2SqVolCu[1]+1043761.529453926*rdx2SqVol[0]*rdx2SqVolSq[1]+1.892833619933881e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-5.4838056e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-1.6744728e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-6.7183704e+8*rdx2SqVolCu[1])-2.72420232e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+9.3076104e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+2979504.0*rdx2SqVolCu[0])*rho[1]+1.032818862000306e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+3.895463326200199e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.565397867170925e+8*rdx2SqVolCu[0]*rho[0])*volFac+((-1.263458491192659e+7*rdx2SqVol[0]*rdx2SqVolCu[1])+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2853825.637446513*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-3.312597052538735e+8*rdx2SqVolR4[1])-1.267455527520777e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1378128.741702675*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(7.845903330673281e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.002634505450664e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.280780073139847e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(1.5631245e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-3.615696e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-967725.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(5.166636929999999e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.9207188e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.240533299999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-3.54553416e+8*phiLy[1]*rdx2SqVolR4[1]+(2.93928705e+8*rdx2SqVol[0]*phiUx[1]-1.35473751e+8*rdx2SqVol[0]*phiLy[1]+((-3.636424649020887e+8*phiUx[0])+5.531752422117667e+8*phiLy[0]-1.163655887686684e+9*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-5.67623952e+8*rdx2SqVolSq[0]*phiUx[1])+3.1293288e+7*rdx2SqVolSq[0]*phiLy[1]+(5.441679977753879e+8*phiUx[0]+2.05578101083412e+8*phiLy[0]-1.799755648262666e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.99762793e+8*rdx2SqVolCu[0]*phiUx[1])+1472823.0*rdx2SqVolCu[0]*phiLy[1]+(3.112358382584934e+8*phiUx[0]-7.738046275219913e+7*phiLy[0]-3.396327436986036e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1.02792888e+8*rdx2SqVolR4[0]*phiUx[1]+(2.064260923741879e+8*bcVals[0]-1.03213046187094e+8*phiUx[0])*rdx2SqVolR4[0]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((6.255090874156799e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.063183084441229e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-1.840632657003176e+8*rdx2SqVolCu[0])*rho[3]+((-1.55944656e+8*rdx2SqVolCu[1])-3.0710148e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.40569768e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.59900048e+8*rdx2SqVolCu[0])*rho[2]+(8.723752800000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+2.6637864e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-3.907165754276457e+7*rho[0]*rdx2SqVolCu[1]-7.578527340329195e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-4.831866111218099e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(1.03700668718898e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.729876357747648e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.580326154677349e+7*rdx2SqVolR4[0])*phiUx[3]+((-4074838.318343807*rdx2SqVol[0]*rdx2SqVolCu[1])-6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.307267512076119e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-2.03852126310076e+8*rdx2SqVolR4[1])-4.004977296097099e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.358473674432715e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.762440955346286e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.04427225e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-1.7179164e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.85606585e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+7448760.0*rdx2SqVolR4[0])*phiUx[2]+(9.268408800000001e+7*rdx2SqVolR4[1]+1.83058407e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-2.64231072e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.13565375e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+1.01926063155038e+8*phiLy[0]*rdx2SqVolR4[1]+(2.5982019e+7*rdx2SqVol[0]*phiUx[1]-5507397.0*rdx2SqVol[0]*phiLy[1]+((-2.616405639024413e+7*phiUx[0])+2.012954270604647e+8*phiLy[0]+5.023498826926871e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.222956e+7*rdx2SqVolSq[0]*phiUx[1]-6.949008e+7*rdx2SqVolSq[0]*phiLy[1]+((-3530369.879035339*phiUx[0])-2.886623335846444e+8*phiLy[0]+2.485380394840879e+8*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2.0477061e+7*rdx2SqVolCu[0]*phiUx[1]+1.43100837e+8*rdx2SqVolCu[0]*phiLy[1]+((-1.793807945138149e+7*phiUx[0])-1.24315031671747e+8*phiLy[0]+1.082621267116283e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[3] = (((5.1981552e+7*rdx2SqVolCu[1]+8.4591528e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.43736648e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-993168.0*rdx2SqVolCu[0])*rho[3]+((-1.773250060898033e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-3.014008121001614e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+5.217992890569751e+7*rdx2SqVolCu[0])*rho[2]+(1.302388584758819e+7*rdx2SqVolCu[1]+470715.9838713786*rdx2SqVol[0]*rdx2SqVolSq[1]+8536308.482054757*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-2.4730888e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-7551544.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-2.2741929e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.292744900000001e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-3.4264296e+7*rdx2SqVolR4[0])*phiUx[3]+((-3.0894696e+7*rdx2SqVolR4[1])-5.9861487e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+705375.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-3.9779376e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.939936e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.7513456e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(2.813584035008861e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+3.39120652485041e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-9.798283298525652e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.440434872903132e+7*rdx2SqVolR4[0])*phiUx[2]+(1155172.233549179*rdx2SqVol[0]*rdx2SqVolCu[1]+1.791344449274544e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3.705960859779652e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-3.397535438501266e+7*phiLy[1]*rdx2SqVolR4[1]+((-5697950.058319833*rdx2SqVol[0]*phiUx[1])-6.553339111300069e+7*rdx2SqVol[0]*phiLy[1]+(7049385.0*phiUx[0]+1561287.0*phiLy[0]+2.2558032e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.623970144356256e+7*rdx2SqVolSq[0]*phiUx[1]+1.15968374092867e+8*rdx2SqVolSq[0]*phiLy[1]+((-1.630608e+7*phiUx[0])+1.969968e+7*phiLy[0]+3.261216e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1287019.405122937*rdx2SqVolCu[0]*phiUx[1]+772143.0538617824*rdx2SqVolCu[0]*phiLy[1]+((-436425.0*phiUx[0])-4.0567527e+7*phiLy[0]+2.4494448e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(2.04493464e+8*rdx2SqVolR4[1]+4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+8.491586400000001e+7*rdx2SqVolR4[0]); 

}

void MGpoissonJacobi2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((6.447991200000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.9688856e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.887905122726076e+7*rdx2SqVolCu[1])-5.601520208069406e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-3.571379299595986e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.214413387187174e+9*rdx2SqVol[0]*rdx2SqVolSq[1])-4.580379954982652e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.840632657003176e+8*rdx2SqVolCu[0])*rho[1]+6.7183704e+8*rho[0]*rdx2SqVolCu[1]+1.307751256e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.97417736e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.59900048e+8*rdx2SqVolCu[0]*rho[0])*volFac+(1.9204101e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.5135219e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-6.07505661e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.513594099999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(5.178539520000001e+8*rdx2SqVolR4[1]+1.008237552e+9*rdx2SqVol[0]*rdx2SqVolCu[1]+2.03977536e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.30827312e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.933865037539783e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2609403.823634816*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.325858046406458e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(3.312597052538735e+8*rdx2SqVolR4[1]+6.446659452024169e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.629314075353864e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.395957580471022e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+3.54553416e+8*phiLy[0]*rdx2SqVolR4[1]+((-4.467607425939947e+8*rdx2SqVol[0]*phiUx[1])-6.504368232599894e+8*rdx2SqVol[0]*phiLy[1]+(4.49890875e+8*phiUx[0]+6.899945009999999e+8*phiLy[0]-8.6379048e+8*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.383275325638488e+8*rdx2SqVolSq[0]*phiUx[1])-2.41723701273902e+8*rdx2SqVolSq[0]*phiLy[1]+(2.1840576e+8*phiUx[0]+1.74784848e+8*phiLy[0]-3.78482016e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2613649.080164168*rdx2SqVolCu[0]*phiUx[1]+9.098581884049787e+7*rdx2SqVolCu[0]*phiLy[1]+((-4.8756675e+7*phiUx[0])-7.904150099999999e+7*phiLy[0]-1.175739312e+9*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+2.580326154677349e+7*rdx2SqVolR4[0]*phiUx[1]+(5.243927039999999e+8*bcVals[0]-7448760.0*phiUx[0])*rdx2SqVolR4[0])/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((2.887905122726076e+7*rdx2SqVolCu[1]+1043761.529453926*rdx2SqVol[0]*rdx2SqVolSq[1]+1.892833619933881e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-5.4838056e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-1.6744728e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-6.7183704e+8*rdx2SqVolCu[1])-2.72420232e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+9.3076104e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+2979504.0*rdx2SqVolCu[0])*rho[1]+1.032818862000306e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+3.895463326200199e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.565397867170925e+8*rdx2SqVolCu[0]*rho[0])*volFac+((-1.263458491192659e+7*rdx2SqVol[0]*rdx2SqVolCu[1])+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2853825.637446513*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-3.312597052538735e+8*rdx2SqVolR4[1])-1.267455527520777e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1378128.741702675*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(7.845903330673281e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.002634505450664e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.280780073139847e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(1.5631245e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-3.615696e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-967725.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(5.166636929999999e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.9207188e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.240533299999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-3.54553416e+8*phiLy[1]*rdx2SqVolR4[1]+(2.93928705e+8*rdx2SqVol[0]*phiUx[1]-1.35473751e+8*rdx2SqVol[0]*phiLy[1]+((-3.636424649020887e+8*phiUx[0])+5.531752422117667e+8*phiLy[0]-1.163655887686684e+9*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-5.67623952e+8*rdx2SqVolSq[0]*phiUx[1])+3.1293288e+7*rdx2SqVolSq[0]*phiLy[1]+(5.441679977753879e+8*phiUx[0]+2.05578101083412e+8*phiLy[0]-1.799755648262666e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.99762793e+8*rdx2SqVolCu[0]*phiUx[1])+1472823.0*rdx2SqVolCu[0]*phiLy[1]+(3.112358382584934e+8*phiUx[0]-7.738046275219913e+7*phiLy[0]-3.396327436986036e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1.02792888e+8*rdx2SqVolR4[0]*phiUx[1]+(2.064260923741879e+8*bcVals[0]-1.03213046187094e+8*phiUx[0])*rdx2SqVolR4[0]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((6.255090874156799e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.063183084441229e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-1.840632657003176e+8*rdx2SqVolCu[0])*rho[3]+((-1.55944656e+8*rdx2SqVolCu[1])-3.0710148e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.40569768e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.59900048e+8*rdx2SqVolCu[0])*rho[2]+(8.723752800000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+2.6637864e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-3.907165754276457e+7*rho[0]*rdx2SqVolCu[1]-7.578527340329195e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-4.831866111218099e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(1.03700668718898e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.729876357747648e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.580326154677349e+7*rdx2SqVolR4[0])*phiUx[3]+((-4074838.318343807*rdx2SqVol[0]*rdx2SqVolCu[1])-6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.307267512076119e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-2.03852126310076e+8*rdx2SqVolR4[1])-4.004977296097099e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.358473674432715e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.762440955346286e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.04427225e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-1.7179164e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.85606585e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+7448760.0*rdx2SqVolR4[0])*phiUx[2]+(9.268408800000001e+7*rdx2SqVolR4[1]+1.83058407e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-2.64231072e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.13565375e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+1.01926063155038e+8*phiLy[0]*rdx2SqVolR4[1]+(2.5982019e+7*rdx2SqVol[0]*phiUx[1]-5507397.0*rdx2SqVol[0]*phiLy[1]+((-2.616405639024413e+7*phiUx[0])+2.012954270604647e+8*phiLy[0]+5.023498826926871e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.222956e+7*rdx2SqVolSq[0]*phiUx[1]-6.949008e+7*rdx2SqVolSq[0]*phiLy[1]+((-3530369.879035339*phiUx[0])-2.886623335846444e+8*phiLy[0]+2.485380394840879e+8*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2.0477061e+7*rdx2SqVolCu[0]*phiUx[1]+1.43100837e+8*rdx2SqVolCu[0]*phiLy[1]+((-1.793807945138149e+7*phiUx[0])-1.24315031671747e+8*phiLy[0]+1.082621267116283e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[3] = (((5.1981552e+7*rdx2SqVolCu[1]+8.4591528e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.43736648e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-993168.0*rdx2SqVolCu[0])*rho[3]+((-1.773250060898033e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-3.014008121001614e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+5.217992890569751e+7*rdx2SqVolCu[0])*rho[2]+(1.302388584758819e+7*rdx2SqVolCu[1]+470715.9838713786*rdx2SqVol[0]*rdx2SqVolSq[1]+8536308.482054757*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-2.4730888e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-7551544.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+((-2.2741929e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.292744900000001e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-3.4264296e+7*rdx2SqVolR4[0])*phiUx[3]+((-3.0894696e+7*rdx2SqVolR4[1])-5.9861487e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+705375.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-3.9779376e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.939936e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.7513456e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(2.813584035008861e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+3.39120652485041e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-9.798283298525652e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.440434872903132e+7*rdx2SqVolR4[0])*phiUx[2]+(1155172.233549179*rdx2SqVol[0]*rdx2SqVolCu[1]+1.791344449274544e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3.705960859779652e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-3.397535438501266e+7*phiLy[1]*rdx2SqVolR4[1]+((-5697950.058319833*rdx2SqVol[0]*phiUx[1])-6.553339111300069e+7*rdx2SqVol[0]*phiLy[1]+(7049385.0*phiUx[0]+1561287.0*phiLy[0]+2.2558032e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.623970144356256e+7*rdx2SqVolSq[0]*phiUx[1]+1.15968374092867e+8*rdx2SqVolSq[0]*phiLy[1]+((-1.630608e+7*phiUx[0])+1.969968e+7*phiLy[0]+3.261216e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1287019.405122937*rdx2SqVolCu[0]*phiUx[1]+772143.0538617824*rdx2SqVolCu[0]*phiLy[1]+((-436425.0*phiUx[0])-4.0567527e+7*phiLy[0]+2.4494448e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])/(2.04493464e+8*rdx2SqVolR4[1]+4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+8.491586400000001e+7*rdx2SqVolR4[0]); 

}

void MGpoissonJacobi2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-40344.0*phiLy[0]*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiUx[1]+50064.92859277839*rdx2SqVol[0]*phiLy[1]+((-52029.0*phiUx[0])-293355.0*phiLy[0]+92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+((-293355.0*phiUx[0])-52029.0*phiLy[0]+828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+(115456.0*bcVals[0]-40344.0*phiUx[0])*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+121032.0*phiLy[1]*rdx2SqVolCu[1]+((-167649.0*rdx2SqVol[0]*phiUx[1])+221031.0*rdx2SqVol[0]*phiLy[1]+(190246.7286525579*phiUx[0]-190246.7286525579*phiLy[0]+373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-11367.0*rdx2SqVolSq[0]*phiUx[1])+29889.0*rdx2SqVolSq[0]*phiLy[1]+(36430.22463559618*phiUx[0]-36430.22463559618*phiLy[0]+1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiUx[1]+136347.039571822*bcVals[0]*rdx2SqVolCu[0])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVol[0]*phiUx[1]+35055.0*rdx2SqVol[0]*phiLy[1]+((-36430.22463559618*phiUx[0])+36430.22463559618*phiLy[0]+64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(188385.0*rdx2SqVolSq[0]*phiUx[1]-35055.0*rdx2SqVolSq[0]*phiLy[1]+((-190246.7286525579*phiUx[0])+190246.7286525579*phiLy[0]+583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLy[2]+((-39128.75979378851*rdx2SqVol[0]*phiUx[1])-190246.7286525579*rdx2SqVol[0]*phiLy[1]+(44403.0*phiUx[0]-44403.0*phiLy[0]+87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(44657.46597154836*rdx2SqVolSq[0]*phiUx[1]-36430.22463559618*rdx2SqVolSq[0]*phiLy[1]+((-44403.0*phiUx[0])+44403.0*phiLy[0]+168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-624456.0*phiUy[0]*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiLx[1]-4159512.0*rdx2SqVol[0]*bcVals[1]+((-3368931.0*phiUy[0])-173313.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiLx[1]-1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+((-173313.0*phiUy[0])-3368931.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiLx[1]-1862784.0*rdx2SqVolCu[0]*bcVals[1]-624456.0*phiLx[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+624456.0*phiUy[1]*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiLx[1]+5603489.203427449*rdx2SqVol[0]*bcVals[1]+((-847040.3948826761*phiUy[0])-1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiLx[1]+5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+((-38955.5547130316*phiUy[0])-2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiLx[1]+733281.0298923596*rdx2SqVolCu[0]*bcVals[1]-366640.5149461798*phiLx[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiLx[1]-934933.3131127583*rdx2SqVol[0]*bcVals[1]+((-2275410.734360502*phiUy[0])-38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiLx[1]-1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+((-1100685.379244677*phiUy[0])-847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiLx[1]+419832.0*rdx2SqVol[0]*bcVals[1]+((-90933.0*phiUy[0])-82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiLx[1]-73032.0*rdx2SqVolSq[0]*bcVals[1]+((-82467.0*phiUy[0])-90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((1.9688856e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+6.447991200000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1.840632657003176e+8*rdx2SqVolCu[1]-4.580379954982652e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.214413387187174e+9*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-3.571379299595986e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.601520208069406e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-2.887905122726076e+7*rdx2SqVolCu[0])*rho[1]-1.59900048e+8*rho[0]*rdx2SqVolCu[1]+2.97417736e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.307751256e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+6.7183704e+8*rdx2SqVolCu[0]*rho[0])*volFac+(1.5135219e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.9204101e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(8.513594099999999e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.07505661e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+(2.580326154677349e+7*rdx2SqVolR4[1]+2613649.080164168*rdx2SqVol[0]*rdx2SqVolCu[1]-3.383275325638488e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-4.467607425939947e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(9.098581884049787e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.41723701273902e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.504368232599894e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(5.243927039999999e+8*rdx2SqVolR4[1]-1.175739312e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-3.78482016e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-8.6379048e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-7448760.0*phiUy[0]*rdx2SqVolR4[1]+((-1.325858046406458e+7*rdx2SqVol[0]*phiUy[1])-7.395957580471022e+7*rdx2SqVol[0]*phiLx[1]-1.30827312e+8*rdx2SqVol[0]*bcVals[1]+((-4.8756675e+7*phiUy[0])-7.904150099999999e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-2609403.823634816*rdx2SqVolSq[0]*phiUy[1])+1.629314075353864e+8*rdx2SqVolSq[0]*phiLx[1]+2.03977536e+8*rdx2SqVolSq[0]*bcVals[1]+(2.1840576e+8*phiUy[0]+1.74784848e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-1.933865037539783e+7*rdx2SqVolCu[0]*phiUy[1])+6.446659452024169e+8*rdx2SqVolCu[0]*phiLx[1]+1.008237552e+9*rdx2SqVolCu[0]*bcVals[1]+(4.49890875e+8*phiUy[0]+6.899945009999999e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+3.312597052538735e+8*rdx2SqVolR4[0]*phiLx[1]+5.178539520000001e+8*rdx2SqVolR4[0]*bcVals[1]+3.54553416e+8*phiLx[0]*rdx2SqVolR4[0])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[1] = (((1.840632657003176e+8*rdx2SqVolCu[1]-1.063183084441229e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-6.255090874156799e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.6637864e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.723752800000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.59900048e+8*rdx2SqVolCu[1])-3.40569768e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.0710148e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.55944656e+8*rdx2SqVolCu[0])*rho[1]+4.831866111218099e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+7.578527340329195e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+3.907165754276457e+7*rdx2SqVolCu[0]*rho[0])*volFac+(2.580326154677349e+7*rdx2SqVolR4[1]+2.729876357747648e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.03700668718898e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1.307267512076119e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4074838.318343807*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-2.0477061e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-1.222956e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.5982019e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-1.43100837e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.949008e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5507397.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-1.082621267116283e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.485380394840879e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-5.023498826926871e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-7448760.0*phiUy[1]*rdx2SqVolR4[1]+((-2.85606585e+8*rdx2SqVol[0]*phiUy[1])+1.13565375e+8*rdx2SqVol[0]*phiLx[1]-1.762440955346286e+8*rdx2SqVol[0]*bcVals[1]+(1.793807945138149e+7*phiUy[0]+1.24315031671747e+8*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.7179164e+8*rdx2SqVolSq[0]*phiUy[1]+2.64231072e+8*rdx2SqVolSq[0]*phiLx[1]-3.358473674432715e+8*rdx2SqVolSq[0]*bcVals[1]+(3530369.879035339*phiUy[0]+2.886623335846444e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.04427225e+8*rdx2SqVolCu[0]*phiUy[1]-1.83058407e+8*rdx2SqVolCu[0]*phiLx[1]+4.004977296097099e+8*rdx2SqVolCu[0]*bcVals[1]+(2.616405639024413e+7*phiUy[0]-2.012954270604647e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-9.268408800000001e+7*rdx2SqVolR4[0]*phiLx[1]+2.03852126310076e+8*rdx2SqVolR4[0]*bcVals[1]-1.01926063155038e+8*phiLx[0]*rdx2SqVolR4[0])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((1.892833619933881e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1043761.529453926*rdx2SqVolSq[0]*rdx2SqVol[1]+2.887905122726076e+7*rdx2SqVolCu[0])*rho[3]+(2979504.0*rdx2SqVolCu[1]+9.3076104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.72420232e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-6.7183704e+8*rdx2SqVolCu[0])*rho[2]+((-1.6744728e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.4838056e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-1.565397867170925e+8*rho[0]*rdx2SqVolCu[1]+3.895463326200199e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.032818862000306e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(2853825.637446513*rdx2SqVol[0]*rdx2SqVolCu[1]+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.263458491192659e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(1378128.741702675*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.267455527520777e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.312597052538735e+8*rdx2SqVolR4[0])*phiLx[3]+(1.02792888e+8*rdx2SqVolR4[1]-2.99762793e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-5.67623952e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.93928705e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(1472823.0*rdx2SqVol[0]*rdx2SqVolCu[1]+3.1293288e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.35473751e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.54553416e+8*rdx2SqVolR4[0])*phiLx[2]+(2.064260923741879e+8*rdx2SqVolR4[1]-3.396327436986036e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.799755648262666e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.163655887686684e+9*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1.03213046187094e+8*phiUy[0]*rdx2SqVolR4[1]+((-967725.0*rdx2SqVol[0]*phiUy[1])-7.240533299999999e+7*rdx2SqVol[0]*phiLx[1]-1.280780073139847e+8*rdx2SqVol[0]*bcVals[1]+(3.112358382584934e+8*phiUy[0]-7.738046275219913e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.615696e+7*rdx2SqVolSq[0]*phiUy[1])+1.9207188e+8*rdx2SqVolSq[0]*phiLx[1]+3.002634505450664e+8*rdx2SqVolSq[0]*bcVals[1]+(5.441679977753879e+8*phiUy[0]+2.05578101083412e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.5631245e+7*rdx2SqVolCu[0]*phiUy[1]+5.166636929999999e+8*rdx2SqVolCu[0]*phiLx[1]+7.845903330673281e+8*rdx2SqVolCu[0]*bcVals[1]+(5.531752422117667e+8*phiLx[0]-3.636424649020887e+8*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[3] = -(1.0*(((993168.0*rdx2SqVolCu[1]+1.43736648e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-8.4591528e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-5.1981552e+7*rdx2SqVolCu[0])*rho[3]+((-8536308.482054757*rdx2SqVol[0]*rdx2SqVolSq[1])-470715.9838713786*rdx2SqVolSq[0]*rdx2SqVol[1]-1.302388584758819e+7*rdx2SqVolCu[0])*rho[2]+((-5.217992890569751e+7*rdx2SqVolCu[1])+3.014008121001614e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.773250060898033e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+7551544.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.4730888e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(3.4264296e+7*rdx2SqVolR4[1]-8.292744900000001e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.2741929e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-705375.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.9861487e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.0894696e+7*rdx2SqVolR4[0])*phiLx[3]+((-1287019.405122937*rdx2SqVol[0]*rdx2SqVolCu[1])-1.623970144356256e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5697950.058319833*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-772143.0538617824*rdx2SqVol[0]*rdx2SqVolCu[1])-1.15968374092867e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+6.553339111300069e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.397535438501266e+7*rdx2SqVolR4[0])*phiLx[2]+((-2.4494448e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.261216e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.2558032e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-3.440434872903132e+7*phiUy[1]*rdx2SqVolR4[1]+(9.798283298525652e+7*rdx2SqVol[0]*phiUy[1]+3.705960859779652e+7*rdx2SqVol[0]*phiLx[1]-5.7513456e+7*rdx2SqVol[0]*bcVals[1]+(436425.0*phiUy[0]+4.0567527e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.39120652485041e+7*rdx2SqVolSq[0]*phiUy[1])-1.791344449274544e+7*rdx2SqVolSq[0]*phiLx[1]+3.939936e+7*rdx2SqVolSq[0]*bcVals[1]+(1.630608e+7*phiUy[0]-1.969968e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.813584035008861e+7*rdx2SqVolCu[0]*phiUy[1])-1155172.233549179*rdx2SqVolCu[0]*phiLx[1]+3.9779376e+7*rdx2SqVolCu[0]*bcVals[1]+((-7049385.0*phiUy[0])-1561287.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(8.491586400000001e+7*rdx2SqVolR4[1]-2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.04493464e+8*rdx2SqVolR4[0]); 

}

void MGpoissonJacobi2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((1.9688856e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+6.447991200000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1.840632657003176e+8*rdx2SqVolCu[1]-4.580379954982652e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.214413387187174e+9*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-3.571379299595986e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.601520208069406e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-2.887905122726076e+7*rdx2SqVolCu[0])*rho[1]-1.59900048e+8*rho[0]*rdx2SqVolCu[1]+2.97417736e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.307751256e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+6.7183704e+8*rdx2SqVolCu[0]*rho[0])*volFac+(1.5135219e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.9204101e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(8.513594099999999e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.07505661e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+(2.580326154677349e+7*rdx2SqVolR4[1]+2613649.080164168*rdx2SqVol[0]*rdx2SqVolCu[1]-3.383275325638488e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-4.467607425939947e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(9.098581884049787e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.41723701273902e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.504368232599894e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(5.243927039999999e+8*rdx2SqVolR4[1]-1.175739312e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-3.78482016e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-8.6379048e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-7448760.0*phiUy[0]*rdx2SqVolR4[1]+((-1.325858046406458e+7*rdx2SqVol[0]*phiUy[1])-7.395957580471022e+7*rdx2SqVol[0]*phiLx[1]-1.30827312e+8*rdx2SqVol[0]*bcVals[1]+((-4.8756675e+7*phiUy[0])-7.904150099999999e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-2609403.823634816*rdx2SqVolSq[0]*phiUy[1])+1.629314075353864e+8*rdx2SqVolSq[0]*phiLx[1]+2.03977536e+8*rdx2SqVolSq[0]*bcVals[1]+(2.1840576e+8*phiUy[0]+1.74784848e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-1.933865037539783e+7*rdx2SqVolCu[0]*phiUy[1])+6.446659452024169e+8*rdx2SqVolCu[0]*phiLx[1]+1.008237552e+9*rdx2SqVolCu[0]*bcVals[1]+(4.49890875e+8*phiUy[0]+6.899945009999999e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+3.312597052538735e+8*rdx2SqVolR4[0]*phiLx[1]+5.178539520000001e+8*rdx2SqVolR4[0]*bcVals[1]+3.54553416e+8*phiLx[0]*rdx2SqVolR4[0])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[1] = (((1.840632657003176e+8*rdx2SqVolCu[1]-1.063183084441229e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-6.255090874156799e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.6637864e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.723752800000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.59900048e+8*rdx2SqVolCu[1])-3.40569768e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.0710148e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.55944656e+8*rdx2SqVolCu[0])*rho[1]+4.831866111218099e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+7.578527340329195e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+3.907165754276457e+7*rdx2SqVolCu[0]*rho[0])*volFac+(2.580326154677349e+7*rdx2SqVolR4[1]+2.729876357747648e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.03700668718898e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1.307267512076119e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4074838.318343807*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-2.0477061e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-1.222956e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.5982019e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-1.43100837e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.949008e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5507397.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-1.082621267116283e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.485380394840879e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-5.023498826926871e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-7448760.0*phiUy[1]*rdx2SqVolR4[1]+((-2.85606585e+8*rdx2SqVol[0]*phiUy[1])+1.13565375e+8*rdx2SqVol[0]*phiLx[1]-1.762440955346286e+8*rdx2SqVol[0]*bcVals[1]+(1.793807945138149e+7*phiUy[0]+1.24315031671747e+8*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.7179164e+8*rdx2SqVolSq[0]*phiUy[1]+2.64231072e+8*rdx2SqVolSq[0]*phiLx[1]-3.358473674432715e+8*rdx2SqVolSq[0]*bcVals[1]+(3530369.879035339*phiUy[0]+2.886623335846444e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.04427225e+8*rdx2SqVolCu[0]*phiUy[1]-1.83058407e+8*rdx2SqVolCu[0]*phiLx[1]+4.004977296097099e+8*rdx2SqVolCu[0]*bcVals[1]+(2.616405639024413e+7*phiUy[0]-2.012954270604647e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-9.268408800000001e+7*rdx2SqVolR4[0]*phiLx[1]+2.03852126310076e+8*rdx2SqVolR4[0]*bcVals[1]-1.01926063155038e+8*phiLx[0]*rdx2SqVolR4[0])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((1.892833619933881e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1043761.529453926*rdx2SqVolSq[0]*rdx2SqVol[1]+2.887905122726076e+7*rdx2SqVolCu[0])*rho[3]+(2979504.0*rdx2SqVolCu[1]+9.3076104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.72420232e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-6.7183704e+8*rdx2SqVolCu[0])*rho[2]+((-1.6744728e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.4838056e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-1.565397867170925e+8*rho[0]*rdx2SqVolCu[1]+3.895463326200199e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.032818862000306e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(2853825.637446513*rdx2SqVol[0]*rdx2SqVolCu[1]+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.263458491192659e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(1378128.741702675*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.267455527520777e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.312597052538735e+8*rdx2SqVolR4[0])*phiLx[3]+(1.02792888e+8*rdx2SqVolR4[1]-2.99762793e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-5.67623952e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.93928705e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(1472823.0*rdx2SqVol[0]*rdx2SqVolCu[1]+3.1293288e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.35473751e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.54553416e+8*rdx2SqVolR4[0])*phiLx[2]+(2.064260923741879e+8*rdx2SqVolR4[1]-3.396327436986036e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.799755648262666e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.163655887686684e+9*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1.03213046187094e+8*phiUy[0]*rdx2SqVolR4[1]+((-967725.0*rdx2SqVol[0]*phiUy[1])-7.240533299999999e+7*rdx2SqVol[0]*phiLx[1]-1.280780073139847e+8*rdx2SqVol[0]*bcVals[1]+(3.112358382584934e+8*phiUy[0]-7.738046275219913e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.615696e+7*rdx2SqVolSq[0]*phiUy[1])+1.9207188e+8*rdx2SqVolSq[0]*phiLx[1]+3.002634505450664e+8*rdx2SqVolSq[0]*bcVals[1]+(5.441679977753879e+8*phiUy[0]+2.05578101083412e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.5631245e+7*rdx2SqVolCu[0]*phiUy[1]+5.166636929999999e+8*rdx2SqVolCu[0]*phiLx[1]+7.845903330673281e+8*rdx2SqVolCu[0]*bcVals[1]+(5.531752422117667e+8*phiLx[0]-3.636424649020887e+8*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[3] = -(1.0*(((993168.0*rdx2SqVolCu[1]+1.43736648e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-8.4591528e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-5.1981552e+7*rdx2SqVolCu[0])*rho[3]+((-8536308.482054757*rdx2SqVol[0]*rdx2SqVolSq[1])-470715.9838713786*rdx2SqVolSq[0]*rdx2SqVol[1]-1.302388584758819e+7*rdx2SqVolCu[0])*rho[2]+((-5.217992890569751e+7*rdx2SqVolCu[1])+3.014008121001614e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.773250060898033e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+7551544.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.4730888e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*volFac+(3.4264296e+7*rdx2SqVolR4[1]-8.292744900000001e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.2741929e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-705375.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.9861487e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.0894696e+7*rdx2SqVolR4[0])*phiLx[3]+((-1287019.405122937*rdx2SqVol[0]*rdx2SqVolCu[1])-1.623970144356256e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5697950.058319833*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-772143.0538617824*rdx2SqVol[0]*rdx2SqVolCu[1])-1.15968374092867e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+6.553339111300069e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.397535438501266e+7*rdx2SqVolR4[0])*phiLx[2]+((-2.4494448e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.261216e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.2558032e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-3.440434872903132e+7*phiUy[1]*rdx2SqVolR4[1]+(9.798283298525652e+7*rdx2SqVol[0]*phiUy[1]+3.705960859779652e+7*rdx2SqVol[0]*phiLx[1]-5.7513456e+7*rdx2SqVol[0]*bcVals[1]+(436425.0*phiUy[0]+4.0567527e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.39120652485041e+7*rdx2SqVolSq[0]*phiUy[1])-1.791344449274544e+7*rdx2SqVolSq[0]*phiLx[1]+3.939936e+7*rdx2SqVolSq[0]*bcVals[1]+(1.630608e+7*phiUy[0]-1.969968e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.813584035008861e+7*rdx2SqVolCu[0]*phiUy[1])-1155172.233549179*rdx2SqVolCu[0]*phiLx[1]+3.9779376e+7*rdx2SqVolCu[0]*bcVals[1]+((-7049385.0*phiUy[0])-1561287.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]))/(8.491586400000001e+7*rdx2SqVolR4[1]-2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.04493464e+8*rdx2SqVolR4[0]); 

}

void MGpoissonJacobi2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-40344.0*phiUy[0]*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiLx[1]-92496.0*rdx2SqVol[0]*bcVals[1]+((-293355.0*phiUy[0])-52029.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiLx[1]-828720.0*rdx2SqVolSq[0]*bcVals[1]+((-52029.0*phiUy[0])-293355.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiLx[1]-115456.0*rdx2SqVolCu[0]*bcVals[1]-40344.0*phiLx[0]*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-121032.0*phiUy[1]*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiLx[1]-373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiLx[1]-1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiLx[1]-136347.039571822*rdx2SqVolCu[0]*bcVals[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiLx[1]-64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiLx[1]-583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiLx[1]-87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiLx[1]-168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiUy[0]-44403.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*volFac+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+624456.0*phiLy[0]*rdx2SqVolCu[1]+((-1544603.07302135*rdx2SqVol[0]*phiLy[1])+100062.3072040616*rdx2SqVol[0]*phiLx[1]+4159512.0*rdx2SqVol[0]*bcVals[1]+(3368931.0*phiLy[0]+173313.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-71036.59977082233*rdx2SqVolSq[0]*phiLy[1])+2951378.203030407*rdx2SqVolSq[0]*phiLx[1]+1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+(173313.0*phiLy[0]+3368931.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiLx[1]+1862784.0*rdx2SqVolCu[0]*bcVals[1]+624456.0*phiLx[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*volFac+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]-624456.0*phiLy[1]*rdx2SqVolCu[1]+((-722367.0*rdx2SqVol[0]*phiLy[1])+1097649.0*rdx2SqVol[0]*phiLx[1]-5603489.203427449*rdx2SqVol[0]*bcVals[1]+(847040.3948826761*phiLy[0]+1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-51597.0*rdx2SqVolSq[0]*phiLy[1])+2182239.0*rdx2SqVolSq[0]*phiLx[1]-5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+(38955.5547130316*phiLy[0]+2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiLx[1]-733281.0298923596*rdx2SqVolCu[0]*bcVals[1]+366640.5149461798*phiLx[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+((-497457.0*rdx2SqVol[0]*phiLy[1])+22491.0*rdx2SqVol[0]*phiLx[1]+934933.3131127583*rdx2SqVol[0]*bcVals[1]+(2275410.734360502*phiLy[0]+38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-451143.0*rdx2SqVolSq[0]*phiLy[1])+757809.0*rdx2SqVolSq[0]*phiLx[1]+1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+(1100685.379244677*phiLy[0]+847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+((-474351.5585164657*rdx2SqVol[0]*phiLy[1])+82239.50439417786*rdx2SqVol[0]*phiLx[1]-419832.0*rdx2SqVol[0]*bcVals[1]+(90933.0*phiLy[0]+82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-109228.3200777161*rdx2SqVolSq[0]*phiLy[1])+82946.18112366594*rdx2SqVolSq[0]*phiLx[1]+73032.0*rdx2SqVolSq[0]*bcVals[1]+(82467.0*phiLy[0]+90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1463496.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-6483052.316323845*rdx2SqVolSq[1])-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])-6483052.316323845*rdx2SqVolSq[0])*rho[1]+1.5082056e+8*rho[0]*rdx2SqVolSq[1]+6.76239224e+8*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+1.5082056e+8*rdx2SqVolSq[0]*rho[0])*volFac+((-1.3788513e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-2998647.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-2998647.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1.3788513e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1.16252928e+8*rdx2SqVolCu[1]+5.22905808e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.2339864e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(7.436442362842056e+7*rdx2SqVolCu[1]+3.323615026097367e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+6.975998306284967e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-3204690.633637355*rdx2SqVol[0]*rdx2SqVolSq[1])-1.476291855196238e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+7.959362400000001e+7*phiLy[0]*rdx2SqVolCu[1]+((-1.476291855196238e+7*rdx2SqVol[0]*phiLy[1])+6.975998306284967e+7*rdx2SqVol[0]*phiLx[1]+1.2339864e+8*rdx2SqVol[0]*bcVals[1]+(3.55706679e+8*phiLy[0]+7.4553345e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-3204690.633637355*rdx2SqVolSq[0]*phiLy[1])+3.323615026097367e+8*rdx2SqVolSq[0]*phiLx[1]+5.22905808e+8*rdx2SqVolSq[0]*bcVals[1]+(7.4553345e+7*phiLy[0]+3.55706679e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+7.436442362842056e+7*rdx2SqVolCu[0]*phiLx[1]+1.16252928e+8*rdx2SqVolCu[0]*bcVals[1]+7.959362400000001e+7*phiLx[0]*rdx2SqVolCu[0])/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((6483052.316323845*rdx2SqVolSq[1]+1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-1.5082056e+8*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-3.5007984e+7*rdx2SqVolSq[0])*rho[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-8771188.427967556*rdx2SqVolSq[0]*rho[0])*volFac+((-7.436442362842056e+7*rdx2SqVolCu[1])-8.604493260170916e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.619246322188773e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-4604440.565570912*rdx2SqVol[0]*rdx2SqVolSq[1])-92486.31697175531*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+((-2.832900731640274e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-7176426.895609816*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-1.8655047e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-4056993.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-5040279.0*rdx2SqVol[0]*rdx2SqVolSq[1])-125001.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]-7.959362400000001e+7*phiLy[1]*rdx2SqVolCu[1]+((-9.1983429e+7*rdx2SqVol[0]*phiLy[1])+1.07116875e+8*rdx2SqVol[0]*phiLx[1]-1.662365553838119e+8*rdx2SqVol[0]*bcVals[1]+(1.172561417439388e+8*phiLx[0]-1.997336039383145e+7*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-1.7305083e+7*rdx2SqVolSq[0]*phiLy[1])+1.13325453e+8*rdx2SqVolSq[0]*phiLx[1]-2.331518580374791e+8*rdx2SqVolSq[0]*bcVals[1]+(1.244999003826421e+8*phiLx[0]-4335757.916097598*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2.0806632e+7*rdx2SqVolCu[0]*phiLx[1]-4.576272223287419e+7*rdx2SqVolCu[0]*bcVals[1]+2.28813611164371e+7*phiLx[0]*rdx2SqVolCu[0]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1]+6483052.316323845*rdx2SqVolSq[0])*rho[3]+((-3.5007984e+7*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-1.5082056e+8*rdx2SqVolSq[0])*rho[2]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-8771188.427967556*rho[0]*rdx2SqVolSq[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-92486.31697175531*rdx2SqVol[0]*rdx2SqVolSq[1])-4604440.565570912*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1.619246322188773e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.604493260170916e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.436442362842056e+7*rdx2SqVolCu[0])*phiLx[3]+((-4.576272223287419e+7*rdx2SqVolCu[1])-2.331518580374791e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.662365553838119e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(2.0806632e+7*rdx2SqVolCu[1]+1.13325453e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.07116875e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-1.7305083e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-9.1983429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.959362400000001e+7*rdx2SqVolCu[0])*phiLx[2]+2.28813611164371e+7*phiLy[0]*rdx2SqVolCu[1]+((-125001.0*rdx2SqVol[0]*phiLy[1])-4056993.0*rdx2SqVol[0]*phiLx[1]-7176426.895609816*rdx2SqVol[0]*bcVals[1]+(1.244999003826421e+8*phiLy[0]-4335757.916097598*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5040279.0*rdx2SqVolSq[0]*phiLy[1])-1.8655047e+7*rdx2SqVolSq[0]*phiLx[1]-2.832900731640274e+7*rdx2SqVolSq[0]*bcVals[1]+(1.172561417439388e+8*phiLy[0]-1.997336039383145e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[3] = (((1.1669328e+7*rdx2SqVolSq[1]+5.2828968e+7*rdx2SqVol[0]*rdx2SqVol[1]+1.1669328e+7*rdx2SqVolSq[0])*rho[3]+(640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1]+2923729.475989186*rdx2SqVolSq[0])*rho[2]+(2923729.475989186*rdx2SqVolSq[1]+640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+892952.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-6935544.0*rdx2SqVolCu[1])-3.7224429e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-8287875.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-8287875.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3.7224429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-6935544.0*rdx2SqVolCu[0])*phiLx[3]+(1436304.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3222576.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-41709.51549706613*rdx2SqVol[0]*rdx2SqVolSq[1])-2076512.411924138*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-9072373.0108449*rdx2SqVol[0]*rdx2SqVolSq[1])-4.075563181606774e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7627120.3721457*rdx2SqVolCu[0])*phiLx[2]-7627120.3721457*phiLy[1]*rdx2SqVolCu[1]+((-4.075563181606774e+7*rdx2SqVol[0]*phiLy[1])-2076512.411924138*rdx2SqVol[0]*phiLx[1]+3222576.0*rdx2SqVol[0]*bcVals[1]+((-56373.0*phiLy[0])-2273067.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9072373.0108449*rdx2SqVolSq[0]*phiLy[1])-41709.51549706613*rdx2SqVolSq[0]*phiLx[1]+1436304.0*rdx2SqVolSq[0]*bcVals[1]+((-2273067.0*phiLy[0])-56373.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(4.5906696e+7*rdx2SqVolCu[1]+2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+4.5906696e+7*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1463496.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-6483052.316323845*rdx2SqVolSq[1])-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])-6483052.316323845*rdx2SqVolSq[0])*rho[1]+1.5082056e+8*rho[0]*rdx2SqVolSq[1]+6.76239224e+8*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+1.5082056e+8*rdx2SqVolSq[0]*rho[0])*volFac+((-1.3788513e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-2998647.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-2998647.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1.3788513e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1.16252928e+8*rdx2SqVolCu[1]+5.22905808e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.2339864e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(7.436442362842056e+7*rdx2SqVolCu[1]+3.323615026097367e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+6.975998306284967e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-3204690.633637355*rdx2SqVol[0]*rdx2SqVolSq[1])-1.476291855196238e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+7.959362400000001e+7*phiLy[0]*rdx2SqVolCu[1]+((-1.476291855196238e+7*rdx2SqVol[0]*phiLy[1])+6.975998306284967e+7*rdx2SqVol[0]*phiLx[1]+1.2339864e+8*rdx2SqVol[0]*bcVals[1]+(3.55706679e+8*phiLy[0]+7.4553345e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-3204690.633637355*rdx2SqVolSq[0]*phiLy[1])+3.323615026097367e+8*rdx2SqVolSq[0]*phiLx[1]+5.22905808e+8*rdx2SqVolSq[0]*bcVals[1]+(7.4553345e+7*phiLy[0]+3.55706679e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+7.436442362842056e+7*rdx2SqVolCu[0]*phiLx[1]+1.16252928e+8*rdx2SqVolCu[0]*bcVals[1]+7.959362400000001e+7*phiLx[0]*rdx2SqVolCu[0])/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((6483052.316323845*rdx2SqVolSq[1]+1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-1.5082056e+8*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-3.5007984e+7*rdx2SqVolSq[0])*rho[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-8771188.427967556*rdx2SqVolSq[0]*rho[0])*volFac+((-7.436442362842056e+7*rdx2SqVolCu[1])-8.604493260170916e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.619246322188773e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-4604440.565570912*rdx2SqVol[0]*rdx2SqVolSq[1])-92486.31697175531*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+((-2.832900731640274e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-7176426.895609816*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-1.8655047e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-4056993.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-5040279.0*rdx2SqVol[0]*rdx2SqVolSq[1])-125001.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]-7.959362400000001e+7*phiLy[1]*rdx2SqVolCu[1]+((-9.1983429e+7*rdx2SqVol[0]*phiLy[1])+1.07116875e+8*rdx2SqVol[0]*phiLx[1]-1.662365553838119e+8*rdx2SqVol[0]*bcVals[1]+(1.172561417439388e+8*phiLx[0]-1.997336039383145e+7*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-1.7305083e+7*rdx2SqVolSq[0]*phiLy[1])+1.13325453e+8*rdx2SqVolSq[0]*phiLx[1]-2.331518580374791e+8*rdx2SqVolSq[0]*bcVals[1]+(1.244999003826421e+8*phiLx[0]-4335757.916097598*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2.0806632e+7*rdx2SqVolCu[0]*phiLx[1]-4.576272223287419e+7*rdx2SqVolCu[0]*bcVals[1]+2.28813611164371e+7*phiLx[0]*rdx2SqVolCu[0]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1]+6483052.316323845*rdx2SqVolSq[0])*rho[3]+((-3.5007984e+7*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-1.5082056e+8*rdx2SqVolSq[0])*rho[2]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-8771188.427967556*rho[0]*rdx2SqVolSq[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-92486.31697175531*rdx2SqVol[0]*rdx2SqVolSq[1])-4604440.565570912*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1.619246322188773e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.604493260170916e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.436442362842056e+7*rdx2SqVolCu[0])*phiLx[3]+((-4.576272223287419e+7*rdx2SqVolCu[1])-2.331518580374791e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.662365553838119e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(2.0806632e+7*rdx2SqVolCu[1]+1.13325453e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.07116875e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-1.7305083e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-9.1983429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.959362400000001e+7*rdx2SqVolCu[0])*phiLx[2]+2.28813611164371e+7*phiLy[0]*rdx2SqVolCu[1]+((-125001.0*rdx2SqVol[0]*phiLy[1])-4056993.0*rdx2SqVol[0]*phiLx[1]-7176426.895609816*rdx2SqVol[0]*bcVals[1]+(1.244999003826421e+8*phiLy[0]-4335757.916097598*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5040279.0*rdx2SqVolSq[0]*phiLy[1])-1.8655047e+7*rdx2SqVolSq[0]*phiLx[1]-2.832900731640274e+7*rdx2SqVolSq[0]*bcVals[1]+(1.172561417439388e+8*phiLy[0]-1.997336039383145e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[3] = (((1.1669328e+7*rdx2SqVolSq[1]+5.2828968e+7*rdx2SqVol[0]*rdx2SqVol[1]+1.1669328e+7*rdx2SqVolSq[0])*rho[3]+(640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1]+2923729.475989186*rdx2SqVolSq[0])*rho[2]+(2923729.475989186*rdx2SqVolSq[1]+640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+892952.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+((-6935544.0*rdx2SqVolCu[1])-3.7224429e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-8287875.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-8287875.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3.7224429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-6935544.0*rdx2SqVolCu[0])*phiLx[3]+(1436304.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3222576.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-41709.51549706613*rdx2SqVol[0]*rdx2SqVolSq[1])-2076512.411924138*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-9072373.0108449*rdx2SqVol[0]*rdx2SqVolSq[1])-4.075563181606774e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7627120.3721457*rdx2SqVolCu[0])*phiLx[2]-7627120.3721457*phiLy[1]*rdx2SqVolCu[1]+((-4.075563181606774e+7*rdx2SqVol[0]*phiLy[1])-2076512.411924138*rdx2SqVol[0]*phiLx[1]+3222576.0*rdx2SqVol[0]*bcVals[1]+((-56373.0*phiLy[0])-2273067.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9072373.0108449*rdx2SqVolSq[0]*phiLy[1])-41709.51549706613*rdx2SqVolSq[0]*phiLx[1]+1436304.0*rdx2SqVolSq[0]*bcVals[1]+((-2273067.0*phiLy[0])-56373.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(4.5906696e+7*rdx2SqVolCu[1]+2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+4.5906696e+7*rdx2SqVolCu[0]); 

}

void MGpoissonJacobi2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*volFac+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+40344.0*phiLy[0]*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiLy[1]+50064.92859277839*rdx2SqVol[0]*phiLx[1]+92496.0*rdx2SqVol[0]*bcVals[1]+(293355.0*phiLy[0]+52029.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+288932.0554646022*rdx2SqVolSq[0]*phiLx[1]+828720.0*rdx2SqVolSq[0]*bcVals[1]+(52029.0*phiLy[0]+293355.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiLx[1]+115456.0*rdx2SqVolCu[0]*bcVals[1]+40344.0*phiLx[0]*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*volFac+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+121032.0*phiLy[1]*rdx2SqVolCu[1]+(221031.0*rdx2SqVol[0]*phiLy[1]-167649.0*rdx2SqVol[0]*phiLx[1]+373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLy[0]-190246.7286525579*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(29889.0*rdx2SqVolSq[0]*phiLy[1]-11367.0*rdx2SqVolSq[0]*phiLx[1]+1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLy[0]-36430.22463559618*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiLx[1]+136347.039571822*rdx2SqVolCu[0]*bcVals[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+(35055.0*rdx2SqVol[0]*phiLy[1]+35055.0*rdx2SqVol[0]*phiLx[1]+64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-35055.0*rdx2SqVolSq[0]*phiLy[1])+188385.0*rdx2SqVolSq[0]*phiLx[1]+583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*volFac+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-190246.7286525579*rdx2SqVol[0]*phiLy[1])-39128.75979378851*rdx2SqVol[0]*phiLx[1]+87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLy[0]-44403.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-36430.22463559618*rdx2SqVolSq[0]*phiLy[1])+44657.46597154836*rdx2SqVolSq[0]*phiLx[1]+168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (16.0*rho[0]*omega*volFac+((-8.660254037844386*rdx2SqVol[1]*phiUy[2])+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]-18.0*phiPrevC[0]+9.0*phiLy[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]-18.0*phiPrevC[0]+9.0*phiLx[0])*rdx2SqVol[0])*omega+18.0*phiPrevC[0]*rdx2SqVol[1]+18.0*phiPrevC[0]*rdx2SqVol[0])/(18.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[1] = (16.0*rho[1]*omega*volFac+((-8.660254037844386*rdx2SqVol[1]*phiUy[3])+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]-18.0*phiPrevC[1]+9.0*phiLy[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-46.0*rdx2SqVol[0]*phiPrevC[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])*omega+18.0*phiPrevC[1]*rdx2SqVol[1]+46.0*rdx2SqVol[0]*phiPrevC[1])/(18.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 
  phiC[2] = (16.0*rho[2]*omega*volFac+((-8.660254037844386*rdx2SqVol[0]*phiUx[3])+8.660254037844386*rdx2SqVol[0]*phiLx[3]-7.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]+((-46.0*rdx2SqVol[1])-18.0*rdx2SqVol[0])*phiPrevC[2]-7.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1])*omega+(46.0*rdx2SqVol[1]+18.0*rdx2SqVol[0])*phiPrevC[2])/(46.0*rdx2SqVol[1]+18.0*rdx2SqVol[0]); 
  phiC[3] = (16.0*rho[3]*omega*volFac+((-7.0*rdx2SqVol[1]*phiUy[3])-7.0*rdx2SqVol[0]*phiUx[3]+((-46.0*rdx2SqVol[1])-46.0*rdx2SqVol[0])*phiPrevC[3]-7.0*rdx2SqVol[1]*phiLy[3]-7.0*rdx2SqVol[0]*phiLx[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1])*omega+(46.0*rdx2SqVol[1]+46.0*rdx2SqVol[0])*phiPrevC[3])/(46.0*rdx2SqVol[1]+46.0*rdx2SqVol[0]); 

}

void MGpoissonDampedJacobi2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[0]*rho[1]+288.0*rho[0]*rdx2SqVol[1]+2096.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(162.0*phiUy[0]-324.0*phiPrevC[0]+162.0*phiLy[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]-31.17691453623978*rdx2SqVol[0]*phiUx[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]+(1179.0*phiUy[0]+54.0*phiUx[0]-3060.0*phiPrevC[0]+1179.0*phiLy[0]+1296.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiUx[1]+(1416.0*phiUx[0]-3528.0*phiPrevC[0]+4224.0*bcVals[0])*rdx2SqVolSq[0])*omega+324.0*phiPrevC[0]*rdx2SqVolSq[1]+3060.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*phiPrevC[0]*rdx2SqVolSq[0])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]+471.1178196587346*rdx2SqVol[0]*rho[0])*omega*volFac+(((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]-324.0*phiPrevC[1]+162.0*phiLy[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]-342.0*rdx2SqVol[0]*phiUx[1]-3060.0*rdx2SqVol[0]*phiPrevC[1]+351.0*rdx2SqVol[0]*phiLy[1]+(265.0037735580381*phiUy[0]+342.9460598986376*phiUx[0]+265.0037735580381*phiLy[0]-1745.907214029428*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]-792.0*rdx2SqVolSq[0]*phiUx[1]-3528.0*rdx2SqVolSq[0]*phiPrevC[1]+(831.384387633061*phiUx[0]-1662.768775266122*bcVals[0])*rdx2SqVolSq[0])*omega+324.0*phiPrevC[1]*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]*phiPrevC[1])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = ((859.0972005541631*rdx2SqVol[0]*rho[3]+(736.0*rdx2SqVol[1]+2096.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(138.0*rdx2SqVol[0]*rdx2SqVol[1]+1416.0*rdx2SqVolSq[0])*phiUx[2]+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiPrevC[2]+((-322.0*rdx2SqVolSq[1])-917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(398.3716857408418*phiUy[0]-398.3716857408418*phiLy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiUy[0]-1134.493278957615*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiPrevC[2])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[0]*rho[2])*omega*volFac+(((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiUx[3]+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiPrevC[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+(876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1]+831.384387633061*rdx2SqVolSq[0])*phiUx[2]-206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiUy[0]-255.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiPrevC[3])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[1]-864.0*rho[0]*rdx2SqVol[1]-2256.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-486.0*phiUy[0])+972.0*phiPrevC[0]-486.0*phiLy[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+467.6537180435967*rdx2SqVol[0]*phiUx[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+((-1269.0*phiUy[0])-486.0*phiUx[0]+3024.0*phiPrevC[0]-1269.0*phiLy[0]+864.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiUx[1]+((-984.0*phiUx[0])+984.0*phiPrevC[0]+2816.0*bcVals[0])*rdx2SqVolSq[0])*omega-972.0*phiPrevC[0]*rdx2SqVolSq[1]-3024.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-984.0*phiPrevC[0]*rdx2SqVolSq[0]))/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]-526.5434455009387*rdx2SqVol[0]*rho[0])*omega*volFac+(((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]-972.0*phiPrevC[1]+486.0*phiLy[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]-522.0*rdx2SqVol[0]*phiUx[1]-3024.0*rdx2SqVol[0]*phiPrevC[1]+243.0*rdx2SqVol[0]*phiLy[1]+((-296.1806880942779*phiUy[0])+592.3613761885558*phiUx[0]-296.1806880942779*phiLy[0]+1163.938142686285*bcVals[0])*rdx2SqVol[0])*rdx2SqVol[1]+24.0*rdx2SqVolSq[0]*phiUx[1]-984.0*rdx2SqVolSq[0]*phiPrevC[1]+1108.512516844081*bcVals[0]*rdx2SqVolSq[0])*omega+972.0*phiPrevC[1]*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]*phiPrevC[1])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((415.6921938165305*rdx2SqVol[0]*rho[3]+((-2208.0*rdx2SqVol[1])-2256.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-1242.0*rdx2SqVol[0]*rdx2SqVol[1])-984.0*rdx2SqVolSq[0])*phiUx[2]+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiPrevC[2]+(966.0*rdx2SqVolSq[1]+987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1195.115057222525*phiLy[0]-1195.115057222525*phiUy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiLy[0]-1221.095819336058*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiPrevC[2]))/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[0]*rho[2])*omega*volFac+(((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUx[3]+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiPrevC[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiLy[0]-285.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiPrevC[3])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[1]-288.0*rho[0]*rdx2SqVol[1]-2096.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(155.8845726811989*rdx2SqVolSq[1]+1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-155.8845726811989*rdx2SqVolSq[1])-1134.493278957615*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-162.0*phiUy[0])+324.0*phiPrevC[0]-162.0*phiLy[0])*rdx2SqVolSq[1]+(483.2421753117166*rdx2SqVol[0]*phiUy[1]+483.2421753117166*rdx2SqVol[0]*phiLy[1]-31.17691453623978*rdx2SqVol[0]*phiLx[1]-1296.0*rdx2SqVol[0]*bcVals[1]+((-1179.0*phiUy[0])+3060.0*phiPrevC[0]-1179.0*phiLy[0]-54.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1247.076581449591*rdx2SqVolSq[0]*phiLx[1]-4224.0*rdx2SqVolSq[0]*bcVals[1]+(3528.0*phiPrevC[0]-1416.0*phiLx[0])*rdx2SqVolSq[0])*omega-324.0*phiPrevC[0]*rdx2SqVolSq[1]-3060.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*phiPrevC[0]*rdx2SqVolSq[0]))/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[1] = (((288.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[1]-471.1178196587346*rdx2SqVol[0]*rho[0])*omega*volFac+(((-155.8845726811989*rdx2SqVolSq[1])-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(155.8845726811989*rdx2SqVolSq[1]+337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-255.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(162.0*phiUy[1]-324.0*phiPrevC[1]+162.0*phiLy[1])*rdx2SqVolSq[1]+(351.0*rdx2SqVol[0]*phiUy[1]-3060.0*rdx2SqVol[0]*phiPrevC[1]+351.0*rdx2SqVol[0]*phiLy[1]-342.0*rdx2SqVol[0]*phiLx[1]+1745.907214029428*rdx2SqVol[0]*bcVals[1]+((-265.0037735580381*phiUy[0])-265.0037735580381*phiLy[0]-342.9460598986376*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0]*phiPrevC[1]-792.0*rdx2SqVolSq[0]*phiLx[1]+1662.768775266122*rdx2SqVolSq[0]*bcVals[1]-831.384387633061*phiLx[0]*rdx2SqVolSq[0])*omega+324.0*phiPrevC[1]*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]*phiPrevC[1])/(324.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[2] = -(1.0*((859.0972005541631*rdx2SqVol[0]*rho[3]+((-736.0*rdx2SqVol[1])-2096.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])-1247.076581449591*rdx2SqVolSq[0])*phiLx[3]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiPrevC[2]+(322.0*rdx2SqVolSq[1]+917.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+((-138.0*rdx2SqVol[0]*rdx2SqVol[1])-1416.0*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiLy[0]-398.3716857408418*phiUy[0])*rdx2SqVolSq[1]+(465.0*rdx2SqVol[0]*phiUy[1]-465.0*rdx2SqVol[0]*phiLy[1]+(1134.493278957615*phiLy[0]-1134.493278957615*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiPrevC[2]))/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 
  phiC[3] = (((736.0*rdx2SqVol[1]+624.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[0]*rho[2])*omega*volFac+(((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-2116.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-3528.0*rdx2SqVolSq[0])*phiPrevC[3]+((-322.0*rdx2SqVolSq[1])-273.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-874.0*rdx2SqVol[0]*rdx2SqVol[1])-792.0*rdx2SqVolSq[0])*phiLx[3]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+206.1140461006964*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+((-876.4177086298519*rdx2SqVol[0]*rdx2SqVol[1])-831.384387633061*rdx2SqVolSq[0])*phiLx[2]+(398.3716857408418*phiUy[1]-398.3716857408418*phiLy[1])*rdx2SqVolSq[1]+(337.749907475931*rdx2SqVol[0]*phiUy[1]-337.749907475931*rdx2SqVol[0]*phiLy[1]+(255.0*phiLy[0]-255.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0])*phiPrevC[3])/(2116.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+3528.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[0]*rho[1]+864.0*rho[0]*rdx2SqVol[1]+2256.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+((-467.6537180435967*rdx2SqVolSq[1])-1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(467.6537180435967*rdx2SqVolSq[1]+1221.095819336058*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(486.0*phiUy[0]-972.0*phiPrevC[0]+486.0*phiLy[0])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]+233.8268590217983*rdx2SqVol[0]*phiLy[1]+467.6537180435967*rdx2SqVol[0]*phiLx[1]+864.0*rdx2SqVol[0]*bcVals[1]+(1269.0*phiUy[0]-3024.0*phiPrevC[0]+1269.0*phiLy[0]+486.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0]*phiLx[1]+2816.0*rdx2SqVolSq[0]*bcVals[1]+(984.0*phiLx[0]-984.0*phiPrevC[0])*rdx2SqVolSq[0])*omega+972.0*phiPrevC[0]*rdx2SqVolSq[1]+3024.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+984.0*phiPrevC[0]*rdx2SqVolSq[0])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[1] = (((864.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[1]+526.5434455009387*rdx2SqVol[0]*rho[0])*omega*volFac+(((-467.6537180435967*rdx2SqVolSq[1])-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+(467.6537180435967*rdx2SqVolSq[1]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]+285.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]+(486.0*phiUy[1]-972.0*phiPrevC[1]+486.0*phiLy[1])*rdx2SqVolSq[1]+(243.0*rdx2SqVol[0]*phiUy[1]-3024.0*rdx2SqVol[0]*phiPrevC[1]+243.0*rdx2SqVol[0]*phiLy[1]-522.0*rdx2SqVol[0]*phiLx[1]+1163.938142686285*rdx2SqVol[0]*bcVals[1]+(296.1806880942779*phiUy[0]+296.1806880942779*phiLy[0]-592.3613761885558*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-984.0*rdx2SqVolSq[0]*phiPrevC[1]+24.0*rdx2SqVolSq[0]*phiLx[1]+1108.512516844081*rdx2SqVolSq[0]*bcVals[1])*omega+972.0*phiPrevC[1]*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]*phiPrevC[1])/(972.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[2] = ((415.6921938165305*rdx2SqVol[0]*rho[3]+(2208.0*rdx2SqVol[1]+2256.0*rdx2SqVol[0])*rho[2])*omega*volFac+((-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[3])-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[3]+(1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1]+969.9484522385712*rdx2SqVolSq[0])*phiLx[3]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiPrevC[2]+((-966.0*rdx2SqVolSq[1])-987.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(1242.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiLx[2]+(1195.115057222525*phiUy[0]-1195.115057222525*phiLy[0])*rdx2SqVolSq[1]+(225.0*rdx2SqVol[0]*phiUy[1]-225.0*rdx2SqVol[0]*phiLy[1]+(1221.095819336058*phiUy[0]-1221.095819336058*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiPrevC[2])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 
  phiC[3] = (((2208.0*rdx2SqVol[1]+432.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[0]*rho[2])*omega*volFac+(((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-6348.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-984.0*rdx2SqVolSq[0])*phiPrevC[3]+((-966.0*rdx2SqVolSq[1])-189.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+(24.0*rdx2SqVolSq[0]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLx[3]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiUy[2]-230.3627574066607*rdx2SqVol[0]*rdx2SqVol[1]*phiLy[2]-1513.812405815199*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(1195.115057222525*phiUy[1]-1195.115057222525*phiLy[1])*rdx2SqVolSq[1]+(233.8268590217983*rdx2SqVol[0]*phiUy[1]-233.8268590217983*rdx2SqVol[0]*phiLy[1]+(285.0*phiUy[0]-285.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0])*phiPrevC[3])/(6348.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+984.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((859.0972005541631*rdx2SqVol[1]*rho[2]+2096.0*rho[0]*rdx2SqVol[1]+288.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(4224.0*rdx2SqVolSq[1]+1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(1416.0*phiUy[0]-3528.0*phiPrevC[0])*rdx2SqVolSq[1]+((-1134.493278957615*rdx2SqVol[0]*phiUx[1])+1134.493278957615*rdx2SqVol[0]*phiLx[1]+(54.0*phiUy[0]+1179.0*phiUx[0]-3060.0*phiPrevC[0]+1179.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-155.8845726811989*rdx2SqVolSq[0]*phiUx[1]+155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+(162.0*phiUx[0]-324.0*phiPrevC[0]+162.0*phiLx[0])*rdx2SqVolSq[0])*omega+3528.0*phiPrevC[0]*rdx2SqVolSq[1]+3060.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+324.0*phiPrevC[0]*rdx2SqVolSq[0])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = ((859.0972005541631*rdx2SqVol[1]*rho[3]+(2096.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[1])*omega*volFac+(((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(1416.0*phiUy[1]-3528.0*phiPrevC[1])*rdx2SqVolSq[1]+(138.0*rdx2SqVol[0]*phiUy[1]-917.0*rdx2SqVol[0]*phiUx[1]-7820.0*rdx2SqVol[0]*phiPrevC[1]-917.0*rdx2SqVol[0]*phiLx[1]+(1134.493278957615*phiUx[0]-1134.493278957615*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-322.0*rdx2SqVolSq[0]*phiUx[1]-2116.0*rdx2SqVolSq[0]*phiPrevC[1]-322.0*rdx2SqVolSq[0]*phiLx[1]+(398.3716857408418*phiUx[0]-398.3716857408418*phiLx[0])*rdx2SqVolSq[0])*omega+3528.0*phiPrevC[1]*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]*phiPrevC[1])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]+471.1178196587346*rho[0]*rdx2SqVol[1])*omega*volFac+(((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+((-3528.0*rdx2SqVolSq[1])-3060.0*rdx2SqVol[0]*rdx2SqVol[1]-324.0*rdx2SqVolSq[0])*phiPrevC[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]+((-1662.768775266122*rdx2SqVolSq[1])-1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+831.384387633061*phiUy[0]*rdx2SqVolSq[1]+((-255.0*rdx2SqVol[0]*phiUx[1])+255.0*rdx2SqVol[0]*phiLx[1]+(342.9460598986376*phiUy[0]+265.0037735580381*phiUx[0]+265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0])*phiPrevC[2])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]+471.1178196587346*rdx2SqVol[1]*rho[1])*omega*volFac+(((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-3528.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-2116.0*rdx2SqVolSq[0])*phiPrevC[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]+831.384387633061*phiUy[1]*rdx2SqVolSq[1]+(876.4177086298519*rdx2SqVol[0]*phiUy[1]-206.1140461006964*rdx2SqVol[0]*phiUx[1]-206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiUx[0]-255.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0])*phiPrevC[3])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[2]-2256.0*rho[0]*rdx2SqVol[1]-864.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(984.0*phiPrevC[0]-984.0*phiUy[0])*rdx2SqVolSq[1]+(1221.095819336058*rdx2SqVol[0]*phiUx[1]-1221.095819336058*rdx2SqVol[0]*phiLx[1]+((-486.0*phiUy[0])-1269.0*phiUx[0]+3024.0*phiPrevC[0]-1269.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0]*phiUx[1]-467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+((-486.0*phiUx[0])+972.0*phiPrevC[0]-486.0*phiLx[0])*rdx2SqVolSq[0])*omega-984.0*phiPrevC[0]*rdx2SqVolSq[1]-3024.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-972.0*phiPrevC[0]*rdx2SqVolSq[0]))/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((415.6921938165305*rdx2SqVol[1]*rho[3]+((-2256.0*rdx2SqVol[1])-2208.0*rdx2SqVol[0])*rho[1])*omega*volFac+((969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(984.0*phiPrevC[1]-984.0*phiUy[1])*rdx2SqVolSq[1]+((-1242.0*rdx2SqVol[0]*phiUy[1])+987.0*rdx2SqVol[0]*phiUx[1]+7728.0*rdx2SqVol[0]*phiPrevC[1]+987.0*rdx2SqVol[0]*phiLx[1]+(1221.095819336058*phiLx[0]-1221.095819336058*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+966.0*rdx2SqVolSq[0]*phiUx[1]+6348.0*rdx2SqVolSq[0]*phiPrevC[1]+966.0*rdx2SqVolSq[0]*phiLx[1]+(1195.115057222525*phiLx[0]-1195.115057222525*phiUx[0])*rdx2SqVolSq[0])*omega-984.0*phiPrevC[1]*rdx2SqVolSq[1]-7728.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]-6348.0*rdx2SqVolSq[0]*phiPrevC[1]))/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]-526.5434455009387*rho[0]*rdx2SqVol[1])*omega*volFac+(((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+((-984.0*rdx2SqVolSq[1])-3024.0*rdx2SqVol[0]*rdx2SqVol[1]-972.0*rdx2SqVolSq[0])*phiPrevC[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[2]+(285.0*rdx2SqVol[0]*phiUx[1]-285.0*rdx2SqVol[0]*phiLx[1]+(592.3613761885558*phiUy[0]-296.1806880942779*phiUx[0]-296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0])*phiPrevC[2])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]-526.5434455009387*rdx2SqVol[1]*rho[1])*omega*volFac+((24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiUy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+((-984.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-6348.0*rdx2SqVolSq[0])*phiPrevC[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+(1513.812405815199*rdx2SqVol[0]*phiUy[1]+230.3627574066607*rdx2SqVol[0]*phiUx[1]+230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiLx[0]-285.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0])*phiPrevC[3])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[2]-2096.0*rho[0]*rdx2SqVol[1]-288.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+((-4224.0*rdx2SqVolSq[1])-1296.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+((-1247.076581449591*rdx2SqVolSq[1])-31.17691453623978*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+483.2421753117166*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(3528.0*phiPrevC[0]-1416.0*phiLy[0])*rdx2SqVolSq[1]+(1134.493278957615*rdx2SqVol[0]*phiUx[1]-1134.493278957615*rdx2SqVol[0]*phiLx[1]+((-1179.0*phiUx[0])+3060.0*phiPrevC[0]-54.0*phiLy[0]-1179.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0]*phiUx[1]-155.8845726811989*rdx2SqVolSq[0]*phiLx[1]+((-162.0*phiUx[0])+324.0*phiPrevC[0]-162.0*phiLx[0])*rdx2SqVolSq[0])*omega-3528.0*phiPrevC[0]*rdx2SqVolSq[1]-3060.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-324.0*phiPrevC[0]*rdx2SqVolSq[0]))/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[1] = -(1.0*((859.0972005541631*rdx2SqVol[1]*rho[3]+((-2096.0*rdx2SqVol[1])-736.0*rdx2SqVol[0])*rho[1])*omega*volFac+((-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+((-1247.076581449591*rdx2SqVolSq[1])-79.67433714816835*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-375.8550252424463*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-465.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(3528.0*phiPrevC[1]-1416.0*phiLy[1])*rdx2SqVolSq[1]+(917.0*rdx2SqVol[0]*phiUx[1]+7820.0*rdx2SqVol[0]*phiPrevC[1]-138.0*rdx2SqVol[0]*phiLy[1]+917.0*rdx2SqVol[0]*phiLx[1]+(1134.493278957615*phiLx[0]-1134.493278957615*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1]+322.0*rdx2SqVolSq[0]*phiUx[1]+2116.0*rdx2SqVolSq[0]*phiPrevC[1]+322.0*rdx2SqVolSq[0]*phiLx[1]+(398.3716857408418*phiLx[0]-398.3716857408418*phiUx[0])*rdx2SqVolSq[0])*omega-3528.0*phiPrevC[1]*rdx2SqVolSq[1]-7820.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]-2116.0*rdx2SqVolSq[0]*phiPrevC[1]))/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 
  phiC[2] = (((624.0*rdx2SqVol[1]+288.0*rdx2SqVol[0])*rho[2]-471.1178196587346*rho[0]*rdx2SqVol[1])*omega*volFac+(((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-155.8845726811989*rdx2SqVolSq[0])*phiUx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+155.8845726811989*rdx2SqVolSq[0])*phiLx[3]+(1662.768775266122*rdx2SqVolSq[1]+1745.907214029428*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiUx[2]+((-3528.0*rdx2SqVolSq[1])-3060.0*rdx2SqVol[0]*rdx2SqVol[1]-324.0*rdx2SqVolSq[0])*phiPrevC[2]+((-792.0*rdx2SqVolSq[1])-342.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(351.0*rdx2SqVol[0]*rdx2SqVol[1]+162.0*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[0]*rdx2SqVolSq[1]+(255.0*rdx2SqVol[0]*phiUx[1]-255.0*rdx2SqVol[0]*phiLx[1]+((-265.0037735580381*phiUx[0])-342.9460598986376*phiLy[0]-265.0037735580381*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0])*phiPrevC[2])/(3528.0*rdx2SqVolSq[1]+3060.0*rdx2SqVol[0]*rdx2SqVol[1]+324.0*rdx2SqVolSq[0]); 
  phiC[3] = (((624.0*rdx2SqVol[1]+736.0*rdx2SqVol[0])*rho[3]-471.1178196587346*rdx2SqVol[1]*rho[1])*omega*volFac+(((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiUx[3]+((-3528.0*rdx2SqVolSq[1])-7820.0*rdx2SqVol[0]*rdx2SqVol[1]-2116.0*rdx2SqVolSq[0])*phiPrevC[3]+((-792.0*rdx2SqVolSq[1])-874.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-273.0*rdx2SqVol[0]*rdx2SqVol[1])-322.0*rdx2SqVolSq[0])*phiLx[3]+(337.749907475931*rdx2SqVol[0]*rdx2SqVol[1]+398.3716857408418*rdx2SqVolSq[0])*phiUx[2]+((-337.749907475931*rdx2SqVol[0]*rdx2SqVol[1])-398.3716857408418*rdx2SqVolSq[0])*phiLx[2]-831.384387633061*phiLy[1]*rdx2SqVolSq[1]+(206.1140461006964*rdx2SqVol[0]*phiUx[1]-876.4177086298519*rdx2SqVol[0]*phiLy[1]+206.1140461006964*rdx2SqVol[0]*phiLx[1]+(255.0*phiLx[0]-255.0*phiUx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0])*phiPrevC[3])/(3528.0*rdx2SqVolSq[1]+7820.0*rdx2SqVol[0]*rdx2SqVol[1]+2116.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((415.6921938165305*rdx2SqVol[1]*rho[2]+2256.0*rho[0]*rdx2SqVol[1]+864.0*rdx2SqVol[0]*rho[0])*omega*volFac+((-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+(2816.0*rdx2SqVolSq[1]+864.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]+(969.9484522385712*rdx2SqVolSq[1]+467.6537180435967*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(984.0*phiLy[0]-984.0*phiPrevC[0])*rdx2SqVolSq[1]+((-1221.095819336058*rdx2SqVol[0]*phiUx[1])+1221.095819336058*rdx2SqVol[0]*phiLx[1]+(1269.0*phiUx[0]-3024.0*phiPrevC[0]+486.0*phiLy[0]+1269.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-467.6537180435967*rdx2SqVolSq[0]*phiUx[1]+467.6537180435967*rdx2SqVolSq[0]*phiLx[1]+(486.0*phiUx[0]-972.0*phiPrevC[0]+486.0*phiLx[0])*rdx2SqVolSq[0])*omega+984.0*phiPrevC[0]*rdx2SqVolSq[1]+3024.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+972.0*phiPrevC[0]*rdx2SqVolSq[0])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[1] = ((415.6921938165305*rdx2SqVol[1]*rho[3]+(2256.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[1])*omega*volFac+((-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[3])+(969.9484522385712*rdx2SqVolSq[1]+1195.115057222525*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]-181.8653347947321*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[3]+225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiUx[2]-225.0*rdx2SqVol[0]*rdx2SqVol[1]*phiLx[2]+(984.0*phiLy[1]-984.0*phiPrevC[1])*rdx2SqVolSq[1]+((-987.0*rdx2SqVol[0]*phiUx[1])-7728.0*rdx2SqVol[0]*phiPrevC[1]+1242.0*rdx2SqVol[0]*phiLy[1]-987.0*rdx2SqVol[0]*phiLx[1]+(1221.095819336058*phiUx[0]-1221.095819336058*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-966.0*rdx2SqVolSq[0]*phiUx[1]-6348.0*rdx2SqVolSq[0]*phiPrevC[1]-966.0*rdx2SqVolSq[0]*phiLx[1]+(1195.115057222525*phiUx[0]-1195.115057222525*phiLx[0])*rdx2SqVolSq[0])*omega+984.0*phiPrevC[1]*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]*phiPrevC[1])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 
  phiC[2] = (((432.0*rdx2SqVol[1]+864.0*rdx2SqVol[0])*rho[2]+526.5434455009387*rho[0]*rdx2SqVol[1])*omega*volFac+(((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-467.6537180435967*rdx2SqVolSq[0])*phiUx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+467.6537180435967*rdx2SqVolSq[0])*phiLx[3]+(1108.512516844081*rdx2SqVolSq[1]+1163.938142686285*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[3]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiUx[2]+((-984.0*rdx2SqVolSq[1])-3024.0*rdx2SqVol[0]*rdx2SqVol[1]-972.0*rdx2SqVolSq[0])*phiPrevC[2]+(24.0*rdx2SqVolSq[1]-522.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[2]+(243.0*rdx2SqVol[0]*rdx2SqVol[1]+486.0*rdx2SqVolSq[0])*phiLx[2]+((-285.0*rdx2SqVol[0]*phiUx[1])+285.0*rdx2SqVol[0]*phiLx[1]+(296.1806880942779*phiUx[0]-592.3613761885558*phiLy[0]+296.1806880942779*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0])*phiPrevC[2])/(984.0*rdx2SqVolSq[1]+3024.0*rdx2SqVol[0]*rdx2SqVol[1]+972.0*rdx2SqVolSq[0]); 
  phiC[3] = (((432.0*rdx2SqVol[1]+2208.0*rdx2SqVol[0])*rho[3]+526.5434455009387*rdx2SqVol[1]*rho[1])*omega*volFac+(((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiUx[3]+((-984.0*rdx2SqVolSq[1])-7728.0*rdx2SqVol[0]*rdx2SqVol[1]-6348.0*rdx2SqVolSq[0])*phiPrevC[3]+(24.0*rdx2SqVolSq[1]-1334.0*rdx2SqVol[0]*rdx2SqVol[1])*phiLy[3]+((-189.0*rdx2SqVol[0]*rdx2SqVol[1])-966.0*rdx2SqVolSq[0])*phiLx[3]+(233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1]+1195.115057222525*rdx2SqVolSq[0])*phiUx[2]+((-233.8268590217983*rdx2SqVol[0]*rdx2SqVol[1])-1195.115057222525*rdx2SqVolSq[0])*phiLx[2]+((-230.3627574066607*rdx2SqVol[0]*phiUx[1])-1513.812405815199*rdx2SqVol[0]*phiLy[1]-230.3627574066607*rdx2SqVol[0]*phiLx[1]+(285.0*phiUx[0]-285.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*omega+(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0])*phiPrevC[3])/(984.0*rdx2SqVolSq[1]+7728.0*rdx2SqVol[0]*rdx2SqVol[1]+6348.0*rdx2SqVolSq[0]); 

}

void MGpoissonDampedJacobi2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(624456.0*phiUy[0]-1555848.0*phiPrevC[0])*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiUx[1]+(3368931.0*phiUy[0]+173313.0*phiUx[0]-1.1189052e+7*phiPrevC[0]+4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiUx[1]+(173313.0*phiUy[0]+3368931.0*phiUx[0]-1.1189052e+7*phiPrevC[0]+1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+(624456.0*phiUx[0]-1555848.0*phiPrevC[0]+1862784.0*bcVals[0])*rdx2SqVolCu[0])*omega+1555848.0*phiPrevC[0]*rdx2SqVolCu[1]+1.1189052e+7*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*phiPrevC[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(624456.0*phiUy[1]-1555848.0*phiPrevC[1])*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1097649.0*rdx2SqVol[0]*phiUx[1]-1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]+(847040.3948826761*phiUy[0]+1100685.379244677*phiUx[0]-5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-2182239.0*rdx2SqVolSq[0]*phiUx[1]-1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]+(38955.5547130316*phiUy[0]+2275410.734360502*phiUx[0]-5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-349272.0*rdx2SqVolCu[0]*phiUx[1]-1555848.0*rdx2SqVolCu[0]*phiPrevC[1]+(366640.5149461798*phiUx[0]-733281.0298923596*bcVals[0])*rdx2SqVolCu[0])*omega+1555848.0*phiPrevC[1]*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]*phiPrevC[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiPrevC[2]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiUx[1]+(2275410.734360502*phiUy[0]+38955.5547130316*phiUx[0]+934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiUx[1]+(1100685.379244677*phiUy[0]+847040.3948826761*phiUx[0]+1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiPrevC[2])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiPrevC[3]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiUx[1]+(90933.0*phiUy[0]+82467.0*phiUx[0]-419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiUx[1]+(82467.0*phiUy[0]+90933.0*phiUx[0]+73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiPrevC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1.1265816e+7*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(1.064828809836548e+7*rdx2SqVolSq[1]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(1.064828809836548e+7*rdx2SqVolSq[0]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-9250416.0*rho[0]*rdx2SqVolSq[1]+1.7580136e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-9250416.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((8660259.0*rdx2SqVol[0]*rdx2SqVolSq[1]-7080939.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(8660259.0*rdx2SqVolSq[0]*rdx2SqVol[1]-7080939.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(1492750.66799516*rdx2SqVolCu[1]-2750120.827394134*rdx2SqVol[0]*rdx2SqVolSq[1]+6151376.711030058*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(7130550.065869739*rdx2SqVol[0]*rdx2SqVolSq[1]-7586460.479438022*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.0336768e+7*rdx2SqVolCu[1]-5.7997776e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1893392e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-430920.0*phiUy[0])-1.4737464e+7*phiPrevC[0])*rdx2SqVolCu[1]+((-7586460.479438022*rdx2SqVol[0]*phiUy[1])+6151376.711030058*rdx2SqVol[0]*phiUx[1]+(711555.0*phiUy[0]-6194475.0*phiUx[0]+2.8535112e+7*phiPrevC[0]+1.1893392e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(7130550.065869739*rdx2SqVolSq[0]*phiUy[1]-2750120.827394134*rdx2SqVolSq[0]*phiUx[1]+((-6194475.0*phiUy[0])+711555.0*phiUx[0]+2.8535112e+7*phiPrevC[0]-5.7997776e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1492750.66799516*rdx2SqVolCu[0]*phiUx[1]+((-430920.0*phiUx[0])-1.4737464e+7*phiPrevC[0]+3.0336768e+7*bcVals[0])*rdx2SqVolCu[0])*omega+1.4737464e+7*phiPrevC[0]*rdx2SqVolCu[1]-2.8535112e+7*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*phiPrevC[0]*rdx2SqVolCu[0])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[1] = (((1.064828809836548e+7*rdx2SqVolSq[1]-1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-9250416.0*rdx2SqVolSq[1])+158424.0*rdx2SqVol[0]*rdx2SqVol[1]-172368.0*rdx2SqVolSq[0])*rho[1]-1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+9056020.719170641*rdx2SqVolSq[0]*rho[0])*omega*volFac+((1492750.66799516*rdx2SqVolCu[1]+4633060.973115181*rdx2SqVol[0]*rdx2SqVolSq[1]+114621.9262924855*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-4658626.043034897*rdx2SqVol[0]*rdx2SqVolSq[1])-1632937.664207363*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(7365267.0*rdx2SqVol[0]*rdx2SqVolSq[1]-6022107.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(5763555.0*rdx2SqVol[0]*rdx2SqVolSq[1]+553725.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.894013253264102e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.164345521036225e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-430920.0*phiUy[1])-1.4737464e+7*phiPrevC[1])*rdx2SqVolCu[1]+((-5756175.0*rdx2SqVol[0]*phiUy[1])+4047057.0*rdx2SqVol[0]*phiUx[1]+2.8535112e+7*rdx2SqVol[0]*phiPrevC[1]+((-6452036.482512711*phiUy[0])-5006934.532233768*phiUx[0]-1.602219050314806e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-115425.0*rdx2SqVolSq[0]*phiUy[1])+1.1487735e+7*rdx2SqVolSq[0]*phiUx[1]+2.8535112e+7*rdx2SqVolSq[0]*phiPrevC[1]+(6064299.588730339*phiUy[0]-1.155226793149618e+7*phiUx[0]+2.261939189589393e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-5946696.0*rdx2SqVolCu[0]*phiUx[1]-1.4737464e+7*rdx2SqVolCu[0]*phiPrevC[1]+(5971002.671980642*phiUx[0]-1.194200534396128e+7*bcVals[0])*rdx2SqVolCu[0])*omega+1.4737464e+7*phiPrevC[1]*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]*phiPrevC[1])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1]-1.064828809836548e+7*rdx2SqVolSq[0])*rho[3]+(172368.0*rdx2SqVolSq[1]-158424.0*rdx2SqVol[0]*rdx2SqVol[1]+9250416.0*rdx2SqVolSq[0])*rho[2]-9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-9056020.719170641*rho[0]*rdx2SqVolSq[1]+1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((1632937.664207363*rdx2SqVol[0]*rdx2SqVolSq[1]+4658626.043034897*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-114621.9262924855*rdx2SqVol[0]*rdx2SqVolSq[1])-4633060.973115181*rdx2SqVolSq[0]*rdx2SqVol[1]-1492750.66799516*rdx2SqVolCu[0])*phiUx[3]+(5946696.0*rdx2SqVolCu[1]-1.1487735e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4047057.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(115425.0*rdx2SqVol[0]*rdx2SqVolSq[1]+5756175.0*rdx2SqVolSq[0]*rdx2SqVol[1]+430920.0*rdx2SqVolCu[0])*phiUx[2]+(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0])*phiPrevC[2]+(1.194200534396128e+7*rdx2SqVolCu[1]-2.261939189589393e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.602219050314806e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-5971002.671980642*phiUy[0]*rdx2SqVolCu[1]+((-553725.0*rdx2SqVol[0]*phiUy[1])+6022107.0*rdx2SqVol[0]*phiUx[1]+(1.155226793149618e+7*phiUy[0]-6064299.588730339*phiUx[0]+1.164345521036225e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5763555.0*rdx2SqVolSq[0]*phiUy[1])-7365267.0*rdx2SqVolSq[0]*phiUx[1]+(5006934.532233768*phiUy[0]+6452036.482512711*phiUx[0]-3.894013253264102e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1.4737464e+7*rdx2SqVolCu[1])+2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1.4737464e+7*rdx2SqVolCu[0])*phiPrevC[2]))/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[3] = -(1.0*(((57456.0*rdx2SqVolSq[1]+3025032.0*rdx2SqVol[0]*rdx2SqVol[1]+57456.0*rdx2SqVolSq[0])*rho[3]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[0])*rho[2]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[1])*rho[1]-2716168.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((1982232.0*rdx2SqVolCu[1]-3365199.0*rdx2SqVol[0]*rdx2SqVolSq[1]-25137.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-25137.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3365199.0*rdx2SqVolSq[0]*rdx2SqVol[1]+1982232.0*rdx2SqVolCu[0])*phiUx[3]+(4912488.0*rdx2SqVolCu[1]-9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]-9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]+4912488.0*rdx2SqVolCu[0])*phiPrevC[3]+(462920.0231865111*rdx2SqVol[0]*rdx2SqVolSq[1]+1320669.688212385*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(31098.97224989918*rdx2SqVol[0]*rdx2SqVolSq[1]+3693399.16129776*rdx2SqVolSq[0]*rdx2SqVol[1]-1990334.223993547*rdx2SqVolCu[0])*phiUx[2]+(8810256.0*rdx2SqVol[0]*rdx2SqVolSq[1]-5228496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-1990334.223993547*phiUy[1]*rdx2SqVolCu[1]+(3693399.16129776*rdx2SqVol[0]*phiUy[1]+1320669.688212385*rdx2SqVol[0]*phiUx[1]+((-156975.0*phiUy[0])-1633905.0*phiUx[0]-5228496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(31098.97224989918*rdx2SqVolSq[0]*phiUy[1]+462920.0231865111*rdx2SqVolSq[0]*phiUx[1]+((-1633905.0*phiUy[0])-156975.0*phiUx[0]+8810256.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-4912488.0*rdx2SqVolCu[1])+9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]-4912488.0*rdx2SqVolCu[0])*phiPrevC[3]))/(4912488.0*rdx2SqVolCu[1]-9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]-9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]+4912488.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1.1265816e+7*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(1.064828809836548e+7*rdx2SqVolSq[1]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(1.064828809836548e+7*rdx2SqVolSq[0]-2.043516497629789e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-9250416.0*rho[0]*rdx2SqVolSq[1]+1.7580136e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-9250416.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((8660259.0*rdx2SqVol[0]*rdx2SqVolSq[1]-7080939.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(8660259.0*rdx2SqVolSq[0]*rdx2SqVol[1]-7080939.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(1492750.66799516*rdx2SqVolCu[1]-2750120.827394134*rdx2SqVol[0]*rdx2SqVolSq[1]+6151376.711030058*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(7130550.065869739*rdx2SqVol[0]*rdx2SqVolSq[1]-7586460.479438022*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.0336768e+7*rdx2SqVolCu[1]-5.7997776e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1893392e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-430920.0*phiUy[0])-1.4737464e+7*phiPrevC[0])*rdx2SqVolCu[1]+((-7586460.479438022*rdx2SqVol[0]*phiUy[1])+6151376.711030058*rdx2SqVol[0]*phiUx[1]+(711555.0*phiUy[0]-6194475.0*phiUx[0]+2.8535112e+7*phiPrevC[0]+1.1893392e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(7130550.065869739*rdx2SqVolSq[0]*phiUy[1]-2750120.827394134*rdx2SqVolSq[0]*phiUx[1]+((-6194475.0*phiUy[0])+711555.0*phiUx[0]+2.8535112e+7*phiPrevC[0]-5.7997776e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1492750.66799516*rdx2SqVolCu[0]*phiUx[1]+((-430920.0*phiUx[0])-1.4737464e+7*phiPrevC[0]+3.0336768e+7*bcVals[0])*rdx2SqVolCu[0])*omega+1.4737464e+7*phiPrevC[0]*rdx2SqVolCu[1]-2.8535112e+7*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*phiPrevC[0]*rdx2SqVolCu[0])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[1] = (((1.064828809836548e+7*rdx2SqVolSq[1]-1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-9250416.0*rdx2SqVolSq[1])+158424.0*rdx2SqVol[0]*rdx2SqVol[1]-172368.0*rdx2SqVolSq[0])*rho[1]-1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+9056020.719170641*rdx2SqVolSq[0]*rho[0])*omega*volFac+((1492750.66799516*rdx2SqVolCu[1]+4633060.973115181*rdx2SqVol[0]*rdx2SqVolSq[1]+114621.9262924855*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-4658626.043034897*rdx2SqVol[0]*rdx2SqVolSq[1])-1632937.664207363*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(7365267.0*rdx2SqVol[0]*rdx2SqVolSq[1]-6022107.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(5763555.0*rdx2SqVol[0]*rdx2SqVolSq[1]+553725.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(3.894013253264102e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.164345521036225e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-430920.0*phiUy[1])-1.4737464e+7*phiPrevC[1])*rdx2SqVolCu[1]+((-5756175.0*rdx2SqVol[0]*phiUy[1])+4047057.0*rdx2SqVol[0]*phiUx[1]+2.8535112e+7*rdx2SqVol[0]*phiPrevC[1]+((-6452036.482512711*phiUy[0])-5006934.532233768*phiUx[0]-1.602219050314806e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-115425.0*rdx2SqVolSq[0]*phiUy[1])+1.1487735e+7*rdx2SqVolSq[0]*phiUx[1]+2.8535112e+7*rdx2SqVolSq[0]*phiPrevC[1]+(6064299.588730339*phiUy[0]-1.155226793149618e+7*phiUx[0]+2.261939189589393e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-5946696.0*rdx2SqVolCu[0]*phiUx[1]-1.4737464e+7*rdx2SqVolCu[0]*phiPrevC[1]+(5971002.671980642*phiUx[0]-1.194200534396128e+7*bcVals[0])*rdx2SqVolCu[0])*omega+1.4737464e+7*phiPrevC[1]*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]*phiPrevC[1])/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1.083065226379279e+7*rdx2SqVol[0]*rdx2SqVol[1]-1.064828809836548e+7*rdx2SqVolSq[0])*rho[3]+(172368.0*rdx2SqVolSq[1]-158424.0*rdx2SqVol[0]*rdx2SqVol[1]+9250416.0*rdx2SqVolSq[0])*rho[2]-9581208.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-9056020.719170641*rho[0]*rdx2SqVolSq[1]+1.73794393723655e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((1632937.664207363*rdx2SqVol[0]*rdx2SqVolSq[1]+4658626.043034897*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-114621.9262924855*rdx2SqVol[0]*rdx2SqVolSq[1])-4633060.973115181*rdx2SqVolSq[0]*rdx2SqVol[1]-1492750.66799516*rdx2SqVolCu[0])*phiUx[3]+(5946696.0*rdx2SqVolCu[1]-1.1487735e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4047057.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(115425.0*rdx2SqVol[0]*rdx2SqVolSq[1]+5756175.0*rdx2SqVolSq[0]*rdx2SqVol[1]+430920.0*rdx2SqVolCu[0])*phiUx[2]+(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0])*phiPrevC[2]+(1.194200534396128e+7*rdx2SqVolCu[1]-2.261939189589393e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.602219050314806e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-5971002.671980642*phiUy[0]*rdx2SqVolCu[1]+((-553725.0*rdx2SqVol[0]*phiUy[1])+6022107.0*rdx2SqVol[0]*phiUx[1]+(1.155226793149618e+7*phiUy[0]-6064299.588730339*phiUx[0]+1.164345521036225e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5763555.0*rdx2SqVolSq[0]*phiUy[1])-7365267.0*rdx2SqVolSq[0]*phiUx[1]+(5006934.532233768*phiUy[0]+6452036.482512711*phiUx[0]-3.894013253264102e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1.4737464e+7*rdx2SqVolCu[1])+2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1.4737464e+7*rdx2SqVolCu[0])*phiPrevC[2]))/(1.4737464e+7*rdx2SqVolCu[1]-2.8535112e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.8535112e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1.4737464e+7*rdx2SqVolCu[0]); 
  phiC[3] = -(1.0*(((57456.0*rdx2SqVolSq[1]+3025032.0*rdx2SqVol[0]*rdx2SqVol[1]+57456.0*rdx2SqVolSq[0])*rho[3]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[0])*rho[2]+(3070371.825561197*rdx2SqVol[0]*rdx2SqVol[1]-3018673.57305688*rdx2SqVolSq[1])*rho[1]-2716168.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((1982232.0*rdx2SqVolCu[1]-3365199.0*rdx2SqVol[0]*rdx2SqVolSq[1]-25137.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-25137.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3365199.0*rdx2SqVolSq[0]*rdx2SqVol[1]+1982232.0*rdx2SqVolCu[0])*phiUx[3]+(4912488.0*rdx2SqVolCu[1]-9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]-9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]+4912488.0*rdx2SqVolCu[0])*phiPrevC[3]+(462920.0231865111*rdx2SqVol[0]*rdx2SqVolSq[1]+1320669.688212385*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(31098.97224989918*rdx2SqVol[0]*rdx2SqVolSq[1]+3693399.16129776*rdx2SqVolSq[0]*rdx2SqVol[1]-1990334.223993547*rdx2SqVolCu[0])*phiUx[2]+(8810256.0*rdx2SqVol[0]*rdx2SqVolSq[1]-5228496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-1990334.223993547*phiUy[1]*rdx2SqVolCu[1]+(3693399.16129776*rdx2SqVol[0]*phiUy[1]+1320669.688212385*rdx2SqVol[0]*phiUx[1]+((-156975.0*phiUy[0])-1633905.0*phiUx[0]-5228496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(31098.97224989918*rdx2SqVolSq[0]*phiUy[1]+462920.0231865111*rdx2SqVolSq[0]*phiUx[1]+((-1633905.0*phiUy[0])-156975.0*phiUx[0]+8810256.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-4912488.0*rdx2SqVolCu[1])+9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]-4912488.0*rdx2SqVolCu[0])*phiPrevC[3]))/(4912488.0*rdx2SqVolCu[1]-9511704.0*rdx2SqVol[0]*rdx2SqVolSq[1]-9511704.0*rdx2SqVolSq[0]*rdx2SqVol[1]+4912488.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiUy[0]-40344.0*phiPrevC[0])*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiUx[1]+(293355.0*phiUy[0]+52029.0*phiUx[0]-345384.0*phiPrevC[0]-92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+(52029.0*phiUy[0]+293355.0*phiUx[0]-345384.0*phiPrevC[0]-828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+(40344.0*phiUx[0]-40344.0*phiPrevC[0]-115456.0*bcVals[0])*rdx2SqVolCu[0])*omega+40344.0*phiPrevC[0]*rdx2SqVolCu[1]+345384.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*phiPrevC[0]*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(121032.0*phiPrevC[1]-121032.0*phiUy[1])*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+167649.0*rdx2SqVol[0]*phiUx[1]+1036152.0*rdx2SqVol[0]*phiPrevC[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiUx[0]-373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+11367.0*rdx2SqVolSq[0]*phiUx[1]+1036152.0*rdx2SqVolSq[0]*phiPrevC[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiUx[0]-1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-2952.0*rdx2SqVolCu[0]*phiUx[1]+121032.0*rdx2SqVolCu[0]*phiPrevC[1]-136347.039571822*bcVals[0]*rdx2SqVolCu[0])*omega-121032.0*phiPrevC[1]*rdx2SqVolCu[1]-1036152.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0]*phiPrevC[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[2]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiUx[1]+((-36430.22463559618*phiUy[0])+36430.22463559618*phiUx[0]-64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiUx[1]+((-190246.7286525579*phiUy[0])+190246.7286525579*phiUx[0]-583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[2]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[3]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUy[2]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiUx[1]+(44403.0*phiUy[0]-44403.0*phiUx[0]-87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiUx[1]+((-44403.0*phiUy[0])+44403.0*phiUx[0]-168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(378861.8654443858*rdx2SqVolSq[1]+2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1]+1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-549960.7724192695*rdx2SqVolCu[1])-2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]-100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(1555848.0*phiPrevC[0]-624456.0*phiLy[0])*rdx2SqVolCu[1]+(100062.3072040616*rdx2SqVol[0]*phiUx[1]-1544603.07302135*rdx2SqVol[0]*phiLy[1]+((-173313.0*phiUx[0])+1.1189052e+7*phiPrevC[0]-3368931.0*phiLy[0]-4159512.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2951378.203030407*rdx2SqVolSq[0]*phiUx[1]-71036.59977082233*rdx2SqVolSq[0]*phiLy[1]+((-3368931.0*phiUx[0])+1.1189052e+7*phiPrevC[0]-173313.0*phiLy[0]-1.1134104e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiUx[1]+((-624456.0*phiUx[0])+1555848.0*phiPrevC[0]-1862784.0*bcVals[0])*rdx2SqVolCu[0])*omega-1555848.0*phiPrevC[0]*rdx2SqVolCu[1]-1.1189052e+7*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*phiPrevC[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(451143.0*rdx2SqVol[0]*rdx2SqVolSq[1]+497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-757809.0*rdx2SqVol[0]*rdx2SqVolSq[1])-22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(1555848.0*phiPrevC[1]-624456.0*phiLy[1])*rdx2SqVolCu[1]+(1097649.0*rdx2SqVol[0]*phiUx[1]+1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]-722367.0*rdx2SqVol[0]*phiLy[1]+((-1100685.379244677*phiUx[0])-847040.3948826761*phiLy[0]+5603489.203427449*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(2182239.0*rdx2SqVolSq[0]*phiUx[1]+1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]-51597.0*rdx2SqVolSq[0]*phiLy[1]+((-2275410.734360502*phiUx[0])-38955.5547130316*phiLy[0]+5563665.891259826*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+349272.0*rdx2SqVolCu[0]*phiUx[1]+1555848.0*rdx2SqVolCu[0]*phiPrevC[1]+(733281.0298923596*bcVals[0]-366640.5149461798*phiUx[0])*rdx2SqVolCu[0])*omega-1555848.0*phiPrevC[1]*rdx2SqVolCu[1]-1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0]*phiPrevC[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = (((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+(275184.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiUx[3]+((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(51597.0*rdx2SqVol[0]*rdx2SqVolSq[1]+722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]+624456.0*rdx2SqVolCu[0])*phiUx[2]+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiPrevC[2]+((-349272.0*rdx2SqVolCu[1])-2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+(22491.0*rdx2SqVol[0]*phiUx[1]-497457.0*rdx2SqVol[0]*phiLy[1]+((-38955.5547130316*phiUx[0])-2275410.734360502*phiLy[0]-934933.3131127583*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(757809.0*rdx2SqVolSq[0]*phiUx[1]-451143.0*rdx2SqVolSq[0]*phiLy[1]+((-847040.3948826761*phiUx[0])-1100685.379244677*phiLy[0]-1708037.655172742*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiPrevC[2])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+(60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1]+69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiUx[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiPrevC[3]+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[3]+(109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1]+474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]+122213.5049820599*rdx2SqVolCu[0])*phiUx[2]+((-82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1])-82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+(82239.50439417786*rdx2SqVol[0]*phiUx[1]-474351.5585164657*rdx2SqVol[0]*phiLy[1]+((-82467.0*phiUx[0])-90933.0*phiLy[0]+419832.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(82946.18112366594*rdx2SqVolSq[0]*phiUx[1]-109228.3200777161*rdx2SqVolSq[0]*phiLy[1]+((-90933.0*phiUx[0])-82467.0*phiLy[0]-73032.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiPrevC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((6.447991200000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.9688856e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.887905122726076e+7*rdx2SqVolCu[1])-5.601520208069406e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-3.571379299595986e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.214413387187174e+9*rdx2SqVol[0]*rdx2SqVolSq[1])-4.580379954982652e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.840632657003176e+8*rdx2SqVolCu[0])*rho[1]+6.7183704e+8*rho[0]*rdx2SqVolCu[1]+1.307751256e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.97417736e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.59900048e+8*rdx2SqVolCu[0]*rho[0])*omega*volFac+((1.9204101e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.5135219e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-6.07505661e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.513594099999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(5.178539520000001e+8*rdx2SqVolR4[1]+1.008237552e+9*rdx2SqVol[0]*rdx2SqVolCu[1]+2.03977536e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.30827312e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.933865037539783e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2609403.823634816*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.325858046406458e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(3.312597052538735e+8*rdx2SqVolR4[1]+6.446659452024169e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.629314075353864e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.395957580471022e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(3.54553416e+8*phiLy[0]-6.13480392e+8*phiPrevC[0])*rdx2SqVolR4[1]+((-4.467607425939947e+8*rdx2SqVol[0]*phiUx[1])-6.504368232599894e+8*rdx2SqVol[0]*phiLy[1]+(4.49890875e+8*phiUx[0]-1.212108912e+9*phiPrevC[0]+6.899945009999999e+8*phiLy[0]-8.6379048e+8*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.383275325638488e+8*rdx2SqVolSq[0]*phiUx[1])-2.41723701273902e+8*rdx2SqVolSq[0]*phiLy[1]+(2.1840576e+8*phiUx[0]+1.397230704e+9*phiPrevC[0]+1.74784848e+8*phiLy[0]-3.78482016e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2613649.080164168*rdx2SqVolCu[0]*phiUx[1]+9.098581884049787e+7*rdx2SqVolCu[0]*phiLy[1]+((-4.8756675e+7*phiUx[0])+7.81081488e+8*phiPrevC[0]-7.904150099999999e+7*phiLy[0]-1.175739312e+9*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+2.580326154677349e+7*rdx2SqVolR4[0]*phiUx[1]+((-7448760.0*phiUx[0])-2.54747592e+8*phiPrevC[0]+5.243927039999999e+8*bcVals[0])*rdx2SqVolR4[0])*omega+6.13480392e+8*phiPrevC[0]*rdx2SqVolR4[1]+1.212108912e+9*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*phiPrevC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*phiPrevC[0]*rdx2SqVolR4[0])/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((2.887905122726076e+7*rdx2SqVolCu[1]+1043761.529453926*rdx2SqVol[0]*rdx2SqVolSq[1]+1.892833619933881e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-5.4838056e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-1.6744728e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-6.7183704e+8*rdx2SqVolCu[1])-2.72420232e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+9.3076104e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+2979504.0*rdx2SqVolCu[0])*rho[1]+1.032818862000306e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+3.895463326200199e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.565397867170925e+8*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-1.263458491192659e+7*rdx2SqVol[0]*rdx2SqVolCu[1])+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2853825.637446513*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-3.312597052538735e+8*rdx2SqVolR4[1])-1.267455527520777e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1378128.741702675*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(7.845903330673281e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.002634505450664e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.280780073139847e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(1.5631245e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-3.615696e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-967725.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(5.166636929999999e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.9207188e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.240533299999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(6.13480392e+8*phiPrevC[1]-3.54553416e+8*phiLy[1])*rdx2SqVolR4[1]+(2.93928705e+8*rdx2SqVol[0]*phiUx[1]+1.212108912e+9*rdx2SqVol[0]*phiPrevC[1]-1.35473751e+8*rdx2SqVol[0]*phiLy[1]+((-3.636424649020887e+8*phiUx[0])+5.531752422117667e+8*phiLy[0]-1.163655887686684e+9*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-5.67623952e+8*rdx2SqVolSq[0]*phiUx[1])-1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]+3.1293288e+7*rdx2SqVolSq[0]*phiLy[1]+(5.441679977753879e+8*phiUx[0]+2.05578101083412e+8*phiLy[0]-1.799755648262666e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.99762793e+8*rdx2SqVolCu[0]*phiUx[1])-7.81081488e+8*rdx2SqVolCu[0]*phiPrevC[1]+1472823.0*rdx2SqVolCu[0]*phiLy[1]+(3.112358382584934e+8*phiUx[0]-7.738046275219913e+7*phiLy[0]-3.396327436986036e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1.02792888e+8*rdx2SqVolR4[0]*phiUx[1]+2.54747592e+8*rdx2SqVolR4[0]*phiPrevC[1]+(2.064260923741879e+8*bcVals[0]-1.03213046187094e+8*phiUx[0])*rdx2SqVolR4[0])*omega-6.13480392e+8*phiPrevC[1]*rdx2SqVolR4[1]-1.212108912e+9*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolCu[1]+1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVolSq[1]+7.81081488e+8*rdx2SqVolCu[0]*phiPrevC[1]*rdx2SqVol[1]-2.54747592e+8*rdx2SqVolR4[0]*phiPrevC[1]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((6.255090874156799e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.063183084441229e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-1.840632657003176e+8*rdx2SqVolCu[0])*rho[3]+((-1.55944656e+8*rdx2SqVolCu[1])-3.0710148e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.40569768e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.59900048e+8*rdx2SqVolCu[0])*rho[2]+(8.723752800000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+2.6637864e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-3.907165754276457e+7*rho[0]*rdx2SqVolCu[1]-7.578527340329195e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-4.831866111218099e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((1.03700668718898e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.729876357747648e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.580326154677349e+7*rdx2SqVolR4[0])*phiUx[3]+((-4074838.318343807*rdx2SqVol[0]*rdx2SqVolCu[1])-6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.307267512076119e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-2.03852126310076e+8*rdx2SqVolR4[1])-4.004977296097099e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.358473674432715e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.762440955346286e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.04427225e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-1.7179164e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.85606585e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+7448760.0*rdx2SqVolR4[0])*phiUx[2]+(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0])*phiPrevC[2]+(9.268408800000001e+7*rdx2SqVolR4[1]+1.83058407e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-2.64231072e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.13565375e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+1.01926063155038e+8*phiLy[0]*rdx2SqVolR4[1]+(2.5982019e+7*rdx2SqVol[0]*phiUx[1]-5507397.0*rdx2SqVol[0]*phiLy[1]+((-2.616405639024413e+7*phiUx[0])+2.012954270604647e+8*phiLy[0]+5.023498826926871e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.222956e+7*rdx2SqVolSq[0]*phiUx[1]-6.949008e+7*rdx2SqVolSq[0]*phiLy[1]+((-3530369.879035339*phiUx[0])-2.886623335846444e+8*phiLy[0]+2.485380394840879e+8*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2.0477061e+7*rdx2SqVolCu[0]*phiUx[1]+1.43100837e+8*rdx2SqVolCu[0]*phiLy[1]+((-1.793807945138149e+7*phiUx[0])-1.24315031671747e+8*phiLy[0]+1.082621267116283e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-6.13480392e+8*rdx2SqVolR4[1])-1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]+1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.54747592e+8*rdx2SqVolR4[0])*phiPrevC[2]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[3] = (((5.1981552e+7*rdx2SqVolCu[1]+8.4591528e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.43736648e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-993168.0*rdx2SqVolCu[0])*rho[3]+((-1.773250060898033e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-3.014008121001614e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+5.217992890569751e+7*rdx2SqVolCu[0])*rho[2]+(1.302388584758819e+7*rdx2SqVolCu[1]+470715.9838713786*rdx2SqVol[0]*rdx2SqVolSq[1]+8536308.482054757*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-2.4730888e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-7551544.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-2.2741929e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.292744900000001e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-3.4264296e+7*rdx2SqVolR4[0])*phiUx[3]+((-2.04493464e+8*rdx2SqVolR4[1])-4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-8.491586400000001e+7*rdx2SqVolR4[0])*phiPrevC[3]+((-3.0894696e+7*rdx2SqVolR4[1])-5.9861487e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+705375.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-3.9779376e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.939936e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.7513456e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(2.813584035008861e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+3.39120652485041e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-9.798283298525652e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.440434872903132e+7*rdx2SqVolR4[0])*phiUx[2]+(1155172.233549179*rdx2SqVol[0]*rdx2SqVolCu[1]+1.791344449274544e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3.705960859779652e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-3.397535438501266e+7*phiLy[1]*rdx2SqVolR4[1]+((-5697950.058319833*rdx2SqVol[0]*phiUx[1])-6.553339111300069e+7*rdx2SqVol[0]*phiLy[1]+(7049385.0*phiUx[0]+1561287.0*phiLy[0]+2.2558032e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.623970144356256e+7*rdx2SqVolSq[0]*phiUx[1]+1.15968374092867e+8*rdx2SqVolSq[0]*phiLy[1]+((-1.630608e+7*phiUx[0])+1.969968e+7*phiLy[0]+3.261216e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1287019.405122937*rdx2SqVolCu[0]*phiUx[1]+772143.0538617824*rdx2SqVolCu[0]*phiLy[1]+((-436425.0*phiUx[0])-4.0567527e+7*phiLy[0]+2.4494448e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(2.04493464e+8*rdx2SqVolR4[1]+4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+8.491586400000001e+7*rdx2SqVolR4[0])*phiPrevC[3])/(2.04493464e+8*rdx2SqVolR4[1]+4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+8.491586400000001e+7*rdx2SqVolR4[0]); 

}

void MGpoissonDampedJacobi2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((6.447991200000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.9688856e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.887905122726076e+7*rdx2SqVolCu[1])-5.601520208069406e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-3.571379299595986e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.214413387187174e+9*rdx2SqVol[0]*rdx2SqVolSq[1])-4.580379954982652e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.840632657003176e+8*rdx2SqVolCu[0])*rho[1]+6.7183704e+8*rho[0]*rdx2SqVolCu[1]+1.307751256e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.97417736e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.59900048e+8*rdx2SqVolCu[0]*rho[0])*omega*volFac+((1.9204101e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.5135219e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-6.07505661e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.513594099999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(5.178539520000001e+8*rdx2SqVolR4[1]+1.008237552e+9*rdx2SqVol[0]*rdx2SqVolCu[1]+2.03977536e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.30827312e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.933865037539783e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2609403.823634816*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.325858046406458e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(3.312597052538735e+8*rdx2SqVolR4[1]+6.446659452024169e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.629314075353864e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.395957580471022e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(3.54553416e+8*phiLy[0]-6.13480392e+8*phiPrevC[0])*rdx2SqVolR4[1]+((-4.467607425939947e+8*rdx2SqVol[0]*phiUx[1])-6.504368232599894e+8*rdx2SqVol[0]*phiLy[1]+(4.49890875e+8*phiUx[0]-1.212108912e+9*phiPrevC[0]+6.899945009999999e+8*phiLy[0]-8.6379048e+8*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.383275325638488e+8*rdx2SqVolSq[0]*phiUx[1])-2.41723701273902e+8*rdx2SqVolSq[0]*phiLy[1]+(2.1840576e+8*phiUx[0]+1.397230704e+9*phiPrevC[0]+1.74784848e+8*phiLy[0]-3.78482016e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2613649.080164168*rdx2SqVolCu[0]*phiUx[1]+9.098581884049787e+7*rdx2SqVolCu[0]*phiLy[1]+((-4.8756675e+7*phiUx[0])+7.81081488e+8*phiPrevC[0]-7.904150099999999e+7*phiLy[0]-1.175739312e+9*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+2.580326154677349e+7*rdx2SqVolR4[0]*phiUx[1]+((-7448760.0*phiUx[0])-2.54747592e+8*phiPrevC[0]+5.243927039999999e+8*bcVals[0])*rdx2SqVolR4[0])*omega+6.13480392e+8*phiPrevC[0]*rdx2SqVolR4[1]+1.212108912e+9*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*phiPrevC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*phiPrevC[0]*rdx2SqVolR4[0])/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[1] = -(1.0*(((2.887905122726076e+7*rdx2SqVolCu[1]+1043761.529453926*rdx2SqVol[0]*rdx2SqVolSq[1]+1.892833619933881e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-5.4838056e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-1.6744728e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-6.7183704e+8*rdx2SqVolCu[1])-2.72420232e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+9.3076104e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+2979504.0*rdx2SqVolCu[0])*rho[1]+1.032818862000306e+9*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+3.895463326200199e+8*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]-1.565397867170925e+8*rdx2SqVolCu[0]*rho[0])*omega*volFac+(((-1.263458491192659e+7*rdx2SqVol[0]*rdx2SqVolCu[1])+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2853825.637446513*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[3]+((-3.312597052538735e+8*rdx2SqVolR4[1])-1.267455527520777e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1378128.741702675*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+(7.845903330673281e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.002634505450664e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.280780073139847e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(1.5631245e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-3.615696e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-967725.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUx[2]+(5.166636929999999e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.9207188e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.240533299999999e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+(6.13480392e+8*phiPrevC[1]-3.54553416e+8*phiLy[1])*rdx2SqVolR4[1]+(2.93928705e+8*rdx2SqVol[0]*phiUx[1]+1.212108912e+9*rdx2SqVol[0]*phiPrevC[1]-1.35473751e+8*rdx2SqVol[0]*phiLy[1]+((-3.636424649020887e+8*phiUx[0])+5.531752422117667e+8*phiLy[0]-1.163655887686684e+9*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-5.67623952e+8*rdx2SqVolSq[0]*phiUx[1])-1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]+3.1293288e+7*rdx2SqVolSq[0]*phiLy[1]+(5.441679977753879e+8*phiUx[0]+2.05578101083412e+8*phiLy[0]-1.799755648262666e+9*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.99762793e+8*rdx2SqVolCu[0]*phiUx[1])-7.81081488e+8*rdx2SqVolCu[0]*phiPrevC[1]+1472823.0*rdx2SqVolCu[0]*phiLy[1]+(3.112358382584934e+8*phiUx[0]-7.738046275219913e+7*phiLy[0]-3.396327436986036e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+1.02792888e+8*rdx2SqVolR4[0]*phiUx[1]+2.54747592e+8*rdx2SqVolR4[0]*phiPrevC[1]+(2.064260923741879e+8*bcVals[0]-1.03213046187094e+8*phiUx[0])*rdx2SqVolR4[0])*omega-6.13480392e+8*phiPrevC[1]*rdx2SqVolR4[1]-1.212108912e+9*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolCu[1]+1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVolSq[1]+7.81081488e+8*rdx2SqVolCu[0]*phiPrevC[1]*rdx2SqVol[1]-2.54747592e+8*rdx2SqVolR4[0]*phiPrevC[1]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((6.255090874156799e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.063183084441229e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-1.840632657003176e+8*rdx2SqVolCu[0])*rho[3]+((-1.55944656e+8*rdx2SqVolCu[1])-3.0710148e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.40569768e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.59900048e+8*rdx2SqVolCu[0])*rho[2]+(8.723752800000001e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+2.6637864e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-3.907165754276457e+7*rho[0]*rdx2SqVolCu[1]-7.578527340329195e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-4.831866111218099e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((1.03700668718898e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.729876357747648e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.580326154677349e+7*rdx2SqVolR4[0])*phiUx[3]+((-4074838.318343807*rdx2SqVol[0]*rdx2SqVolCu[1])-6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.307267512076119e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-2.03852126310076e+8*rdx2SqVolR4[1])-4.004977296097099e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+3.358473674432715e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.762440955346286e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+((-1.04427225e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-1.7179164e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.85606585e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+7448760.0*rdx2SqVolR4[0])*phiUx[2]+(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0])*phiPrevC[2]+(9.268408800000001e+7*rdx2SqVolR4[1]+1.83058407e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-2.64231072e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.13565375e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]+1.01926063155038e+8*phiLy[0]*rdx2SqVolR4[1]+(2.5982019e+7*rdx2SqVol[0]*phiUx[1]-5507397.0*rdx2SqVol[0]*phiLy[1]+((-2.616405639024413e+7*phiUx[0])+2.012954270604647e+8*phiLy[0]+5.023498826926871e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.222956e+7*rdx2SqVolSq[0]*phiUx[1]-6.949008e+7*rdx2SqVolSq[0]*phiLy[1]+((-3530369.879035339*phiUx[0])-2.886623335846444e+8*phiLy[0]+2.485380394840879e+8*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(2.0477061e+7*rdx2SqVolCu[0]*phiUx[1]+1.43100837e+8*rdx2SqVolCu[0]*phiLy[1]+((-1.793807945138149e+7*phiUx[0])-1.24315031671747e+8*phiLy[0]+1.082621267116283e+8*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-6.13480392e+8*rdx2SqVolR4[1])-1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]+1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.54747592e+8*rdx2SqVolR4[0])*phiPrevC[2]))/(6.13480392e+8*rdx2SqVolR4[1]+1.212108912e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-7.81081488e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.54747592e+8*rdx2SqVolR4[0]); 
  phiC[3] = (((5.1981552e+7*rdx2SqVolCu[1]+8.4591528e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.43736648e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-993168.0*rdx2SqVolCu[0])*rho[3]+((-1.773250060898033e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-3.014008121001614e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+5.217992890569751e+7*rdx2SqVolCu[0])*rho[2]+(1.302388584758819e+7*rdx2SqVolCu[1]+470715.9838713786*rdx2SqVol[0]*rdx2SqVolSq[1]+8536308.482054757*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-2.4730888e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]-7551544.0*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-2.2741929e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+8.292744900000001e+7*rdx2SqVolCu[0]*rdx2SqVol[1]-3.4264296e+7*rdx2SqVolR4[0])*phiUx[3]+((-2.04493464e+8*rdx2SqVolR4[1])-4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-8.491586400000001e+7*rdx2SqVolR4[0])*phiPrevC[3]+((-3.0894696e+7*rdx2SqVolR4[1])-5.9861487e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+705375.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[3]+((-3.9779376e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.939936e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.7513456e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[3]+(2.813584035008861e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+3.39120652485041e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-9.798283298525652e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.440434872903132e+7*rdx2SqVolR4[0])*phiUx[2]+(1155172.233549179*rdx2SqVol[0]*rdx2SqVolCu[1]+1.791344449274544e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-3.705960859779652e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLy[2]-3.397535438501266e+7*phiLy[1]*rdx2SqVolR4[1]+((-5697950.058319833*rdx2SqVol[0]*phiUx[1])-6.553339111300069e+7*rdx2SqVol[0]*phiLy[1]+(7049385.0*phiUx[0]+1561287.0*phiLy[0]+2.2558032e+7*bcVals[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.623970144356256e+7*rdx2SqVolSq[0]*phiUx[1]+1.15968374092867e+8*rdx2SqVolSq[0]*phiLy[1]+((-1.630608e+7*phiUx[0])+1.969968e+7*phiLy[0]+3.261216e+7*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1287019.405122937*rdx2SqVolCu[0]*phiUx[1]+772143.0538617824*rdx2SqVolCu[0]*phiLy[1]+((-436425.0*phiUx[0])-4.0567527e+7*phiLy[0]+2.4494448e+7*bcVals[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+(2.04493464e+8*rdx2SqVolR4[1]+4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+8.491586400000001e+7*rdx2SqVolR4[0])*phiPrevC[3])/(2.04493464e+8*rdx2SqVolR4[1]+4.04036304e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.60360496e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+8.491586400000001e+7*rdx2SqVolR4[0]); 

}

void MGpoissonDampedJacobi2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-17043.37994647775*rdx2SqVolSq[1])-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[3]+(49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-115456.0*rdx2SqVolCu[1])-828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]-92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1])-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-39767.88654178142*rdx2SqVolCu[1])-288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]-50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(40344.0*phiPrevC[0]-40344.0*phiLy[0])*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiUx[1]+50064.92859277839*rdx2SqVol[0]*phiLy[1]+((-52029.0*phiUx[0])+345384.0*phiPrevC[0]-293355.0*phiLy[0]+92496.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(288932.0554646022*rdx2SqVolSq[0]*phiUx[1]+9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+((-293355.0*phiUx[0])+345384.0*phiPrevC[0]-52029.0*phiLy[0]+828720.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiUx[1]+((-40344.0*phiUx[0])+40344.0*phiPrevC[0]+115456.0*bcVals[0])*rdx2SqVolCu[0])*omega-40344.0*phiPrevC[0]*rdx2SqVolCu[1]-345384.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-345384.0*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-40344.0*phiPrevC[0]*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiUx[3]+(119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1])-64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(35055.0*rdx2SqVol[0]*rdx2SqVolSq[1]-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+((-188385.0*rdx2SqVol[0]*rdx2SqVolSq[1])-35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(121032.0*phiLy[1]-121032.0*phiPrevC[1])*rdx2SqVolCu[1]+((-167649.0*rdx2SqVol[0]*phiUx[1])-1036152.0*rdx2SqVol[0]*phiPrevC[1]+221031.0*rdx2SqVol[0]*phiLy[1]+(190246.7286525579*phiUx[0]-190246.7286525579*phiLy[0]+373818.1334927453*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-11367.0*rdx2SqVolSq[0]*phiUx[1])-1036152.0*rdx2SqVolSq[0]*phiPrevC[1]+29889.0*rdx2SqVolSq[0]*phiLy[1]+(36430.22463559618*phiUx[0]-36430.22463559618*phiLy[0]+1029337.010328493*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0]*phiUx[1]-121032.0*rdx2SqVolCu[0]*phiPrevC[1]+136347.039571822*bcVals[0]*rdx2SqVolCu[0])*omega+121032.0*phiPrevC[1]*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]*phiPrevC[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+((-53136.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiUx[3]+(35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-136347.039571822*rdx2SqVolCu[1])-1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]-373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-29889.0*rdx2SqVol[0]*rdx2SqVolSq[1])-221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiUx[2]+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[2]+((-2952.0*rdx2SqVolCu[1])+11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]+167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVol[0]*phiUx[1]+35055.0*rdx2SqVol[0]*phiLy[1]+((-36430.22463559618*phiUx[0])+36430.22463559618*phiLy[0]+64764.84379661545*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(188385.0*rdx2SqVolSq[0]*phiUx[1]-35055.0*rdx2SqVolSq[0]*phiLy[1]+((-190246.7286525579*phiUx[0])+190246.7286525579*phiLy[0]+583936.681060541*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[2]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+((-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])-64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiUx[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[3]+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-168112.0*rdx2SqVol[0]*rdx2SqVolSq[1])-87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1]+190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUx[2]+(39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1]-44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLy[2]+((-39128.75979378851*rdx2SqVol[0]*phiUx[1])-190246.7286525579*rdx2SqVol[0]*phiLy[1]+(44403.0*phiUx[0]-44403.0*phiLy[0]+87248.0*bcVals[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(44657.46597154836*rdx2SqVolSq[0]*phiUx[1]-36430.22463559618*rdx2SqVolSq[0]*phiLy[1]+((-44403.0*phiUx[0])+44403.0*phiLy[0]+168112.0*bcVals[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[1]-924336.0*rho[0]*rdx2SqVolSq[1]-5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1862784.0*rdx2SqVolCu[1])-1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(1555848.0*phiPrevC[0]-624456.0*phiUy[0])*rdx2SqVolCu[1]+(1544603.07302135*rdx2SqVol[0]*phiUy[1]-100062.3072040616*rdx2SqVol[0]*phiLx[1]-4159512.0*rdx2SqVol[0]*bcVals[1]+((-3368931.0*phiUy[0])+1.1189052e+7*phiPrevC[0]-173313.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(71036.59977082233*rdx2SqVolSq[0]*phiUy[1]-2951378.203030407*rdx2SqVolSq[0]*phiLx[1]-1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+((-173313.0*phiUy[0])+1.1189052e+7*phiPrevC[0]-3368931.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0]*phiLx[1]-1862784.0*rdx2SqVolCu[0]*bcVals[1]+(1555848.0*phiPrevC[0]-624456.0*phiLx[0])*rdx2SqVolCu[0])*omega-1555848.0*phiPrevC[0]*rdx2SqVolCu[1]-1.1189052e+7*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*phiPrevC[0]*rdx2SqVolCu[0]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = (((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(924336.0*rdx2SqVolSq[1]+1737060.0*rdx2SqVol[0]*rdx2SqVol[1]+275184.0*rdx2SqVolSq[0])*rho[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1])-934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(624456.0*phiUy[1]-1555848.0*phiPrevC[1])*rdx2SqVolCu[1]+(722367.0*rdx2SqVol[0]*phiUy[1]-1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]-1097649.0*rdx2SqVol[0]*phiLx[1]+5603489.203427449*rdx2SqVol[0]*bcVals[1]+((-847040.3948826761*phiUy[0])-1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(51597.0*rdx2SqVolSq[0]*phiUy[1]-1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]-2182239.0*rdx2SqVolSq[0]*phiLx[1]+5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+((-38955.5547130316*phiUy[0])-2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0]*phiPrevC[1]-349272.0*rdx2SqVolCu[0]*phiLx[1]+733281.0298923596*rdx2SqVolCu[0]*bcVals[1]-366640.5149461798*phiLx[0]*rdx2SqVolCu[0])*omega+1555848.0*phiPrevC[1]*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]*phiPrevC[1])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]+537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-207762.9584695019*rho[0]*rdx2SqVolSq[1]-1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiPrevC[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+(733281.0298923596*rdx2SqVolCu[1]+5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]+5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]-366640.5149461798*phiUy[0]*rdx2SqVolCu[1]+(497457.0*rdx2SqVol[0]*phiUy[1]-22491.0*rdx2SqVol[0]*phiLx[1]-934933.3131127583*rdx2SqVol[0]*bcVals[1]+((-2275410.734360502*phiUy[0])-38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(451143.0*rdx2SqVolSq[0]*phiUy[1]-757809.0*rdx2SqVolSq[0]*phiLx[1]-1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+((-1100685.379244677*phiUy[0])-847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiPrevC[2]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+(69254.31948983399*rdx2SqVolSq[1]+60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiPrevC[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]+(419832.0*rdx2SqVolSq[0]*rdx2SqVol[1]-73032.0*rdx2SqVol[0]*rdx2SqVolSq[1])*bcVals[2]+122213.5049820599*phiUy[1]*rdx2SqVolCu[1]+(474351.5585164657*rdx2SqVol[0]*phiUy[1]-82239.50439417786*rdx2SqVol[0]*phiLx[1]+419832.0*rdx2SqVol[0]*bcVals[1]+((-90933.0*phiUy[0])-82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(109228.3200777161*rdx2SqVolSq[0]*phiUy[1]-82946.18112366594*rdx2SqVolSq[0]*phiLx[1]-73032.0*rdx2SqVolSq[0]*bcVals[1]+((-82467.0*phiUy[0])-90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiPrevC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((1.9688856e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+6.447991200000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1.840632657003176e+8*rdx2SqVolCu[1]-4.580379954982652e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.214413387187174e+9*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-3.571379299595986e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.601520208069406e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-2.887905122726076e+7*rdx2SqVolCu[0])*rho[1]-1.59900048e+8*rho[0]*rdx2SqVolCu[1]+2.97417736e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.307751256e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+6.7183704e+8*rdx2SqVolCu[0]*rho[0])*omega*volFac+((1.5135219e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.9204101e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(8.513594099999999e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.07505661e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+(2.580326154677349e+7*rdx2SqVolR4[1]+2613649.080164168*rdx2SqVol[0]*rdx2SqVolCu[1]-3.383275325638488e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-4.467607425939947e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(9.098581884049787e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.41723701273902e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.504368232599894e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(5.243927039999999e+8*rdx2SqVolR4[1]-1.175739312e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-3.78482016e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-8.6379048e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-7448760.0*phiUy[0])-2.54747592e+8*phiPrevC[0])*rdx2SqVolR4[1]+((-1.325858046406458e+7*rdx2SqVol[0]*phiUy[1])-7.395957580471022e+7*rdx2SqVol[0]*phiLx[1]-1.30827312e+8*rdx2SqVol[0]*bcVals[1]+((-4.8756675e+7*phiUy[0])+7.81081488e+8*phiPrevC[0]-7.904150099999999e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-2609403.823634816*rdx2SqVolSq[0]*phiUy[1])+1.629314075353864e+8*rdx2SqVolSq[0]*phiLx[1]+2.03977536e+8*rdx2SqVolSq[0]*bcVals[1]+(2.1840576e+8*phiUy[0]+1.397230704e+9*phiPrevC[0]+1.74784848e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-1.933865037539783e+7*rdx2SqVolCu[0]*phiUy[1])+6.446659452024169e+8*rdx2SqVolCu[0]*phiLx[1]+1.008237552e+9*rdx2SqVolCu[0]*bcVals[1]+(4.49890875e+8*phiUy[0]-1.212108912e+9*phiPrevC[0]+6.899945009999999e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+3.312597052538735e+8*rdx2SqVolR4[0]*phiLx[1]+5.178539520000001e+8*rdx2SqVolR4[0]*bcVals[1]+(3.54553416e+8*phiLx[0]-6.13480392e+8*phiPrevC[0])*rdx2SqVolR4[0])*omega+2.54747592e+8*phiPrevC[0]*rdx2SqVolR4[1]-7.81081488e+8*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*phiPrevC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*phiPrevC[0]*rdx2SqVolR4[0])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[1] = (((1.840632657003176e+8*rdx2SqVolCu[1]-1.063183084441229e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-6.255090874156799e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.6637864e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.723752800000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.59900048e+8*rdx2SqVolCu[1])-3.40569768e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.0710148e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.55944656e+8*rdx2SqVolCu[0])*rho[1]+4.831866111218099e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+7.578527340329195e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+3.907165754276457e+7*rdx2SqVolCu[0]*rho[0])*omega*volFac+((2.580326154677349e+7*rdx2SqVolR4[1]+2.729876357747648e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.03700668718898e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1.307267512076119e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4074838.318343807*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-2.0477061e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-1.222956e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.5982019e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-1.43100837e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.949008e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5507397.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-1.082621267116283e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.485380394840879e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-5.023498826926871e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-7448760.0*phiUy[1])-2.54747592e+8*phiPrevC[1])*rdx2SqVolR4[1]+((-2.85606585e+8*rdx2SqVol[0]*phiUy[1])+7.81081488e+8*rdx2SqVol[0]*phiPrevC[1]+1.13565375e+8*rdx2SqVol[0]*phiLx[1]-1.762440955346286e+8*rdx2SqVol[0]*bcVals[1]+(1.793807945138149e+7*phiUy[0]+1.24315031671747e+8*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.7179164e+8*rdx2SqVolSq[0]*phiUy[1]+1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]+2.64231072e+8*rdx2SqVolSq[0]*phiLx[1]-3.358473674432715e+8*rdx2SqVolSq[0]*bcVals[1]+(3530369.879035339*phiUy[0]+2.886623335846444e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.04427225e+8*rdx2SqVolCu[0]*phiUy[1]-1.212108912e+9*rdx2SqVolCu[0]*phiPrevC[1]-1.83058407e+8*rdx2SqVolCu[0]*phiLx[1]+4.004977296097099e+8*rdx2SqVolCu[0]*bcVals[1]+(2.616405639024413e+7*phiUy[0]-2.012954270604647e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-6.13480392e+8*rdx2SqVolR4[0]*phiPrevC[1]-9.268408800000001e+7*rdx2SqVolR4[0]*phiLx[1]+2.03852126310076e+8*rdx2SqVolR4[0]*bcVals[1]-1.01926063155038e+8*phiLx[0]*rdx2SqVolR4[0])*omega+2.54747592e+8*phiPrevC[1]*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*phiPrevC[1]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]*phiPrevC[1])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((1.892833619933881e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1043761.529453926*rdx2SqVolSq[0]*rdx2SqVol[1]+2.887905122726076e+7*rdx2SqVolCu[0])*rho[3]+(2979504.0*rdx2SqVolCu[1]+9.3076104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.72420232e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-6.7183704e+8*rdx2SqVolCu[0])*rho[2]+((-1.6744728e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.4838056e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-1.565397867170925e+8*rho[0]*rdx2SqVolCu[1]+3.895463326200199e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.032818862000306e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2853825.637446513*rdx2SqVol[0]*rdx2SqVolCu[1]+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.263458491192659e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(1378128.741702675*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.267455527520777e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.312597052538735e+8*rdx2SqVolR4[0])*phiLx[3]+(1.02792888e+8*rdx2SqVolR4[1]-2.99762793e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-5.67623952e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.93928705e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0])*phiPrevC[2]+(1472823.0*rdx2SqVol[0]*rdx2SqVolCu[1]+3.1293288e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.35473751e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.54553416e+8*rdx2SqVolR4[0])*phiLx[2]+(2.064260923741879e+8*rdx2SqVolR4[1]-3.396327436986036e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.799755648262666e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.163655887686684e+9*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1.03213046187094e+8*phiUy[0]*rdx2SqVolR4[1]+((-967725.0*rdx2SqVol[0]*phiUy[1])-7.240533299999999e+7*rdx2SqVol[0]*phiLx[1]-1.280780073139847e+8*rdx2SqVol[0]*bcVals[1]+(3.112358382584934e+8*phiUy[0]-7.738046275219913e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.615696e+7*rdx2SqVolSq[0]*phiUy[1])+1.9207188e+8*rdx2SqVolSq[0]*phiLx[1]+3.002634505450664e+8*rdx2SqVolSq[0]*bcVals[1]+(5.441679977753879e+8*phiUy[0]+2.05578101083412e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.5631245e+7*rdx2SqVolCu[0]*phiUy[1]+5.166636929999999e+8*rdx2SqVolCu[0]*phiLx[1]+7.845903330673281e+8*rdx2SqVolCu[0]*bcVals[1]+(5.531752422117667e+8*phiLx[0]-3.636424649020887e+8*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-2.54747592e+8*rdx2SqVolR4[1])+7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]-6.13480392e+8*rdx2SqVolR4[0])*phiPrevC[2]))/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[3] = -(1.0*(((993168.0*rdx2SqVolCu[1]+1.43736648e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-8.4591528e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-5.1981552e+7*rdx2SqVolCu[0])*rho[3]+((-8536308.482054757*rdx2SqVol[0]*rdx2SqVolSq[1])-470715.9838713786*rdx2SqVolSq[0]*rdx2SqVol[1]-1.302388584758819e+7*rdx2SqVolCu[0])*rho[2]+((-5.217992890569751e+7*rdx2SqVolCu[1])+3.014008121001614e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.773250060898033e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+7551544.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.4730888e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((3.4264296e+7*rdx2SqVolR4[1]-8.292744900000001e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.2741929e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(8.491586400000001e+7*rdx2SqVolR4[1]-2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.04493464e+8*rdx2SqVolR4[0])*phiPrevC[3]+((-705375.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.9861487e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.0894696e+7*rdx2SqVolR4[0])*phiLx[3]+((-1287019.405122937*rdx2SqVol[0]*rdx2SqVolCu[1])-1.623970144356256e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5697950.058319833*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-772143.0538617824*rdx2SqVol[0]*rdx2SqVolCu[1])-1.15968374092867e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+6.553339111300069e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.397535438501266e+7*rdx2SqVolR4[0])*phiLx[2]+((-2.4494448e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.261216e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.2558032e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-3.440434872903132e+7*phiUy[1]*rdx2SqVolR4[1]+(9.798283298525652e+7*rdx2SqVol[0]*phiUy[1]+3.705960859779652e+7*rdx2SqVol[0]*phiLx[1]-5.7513456e+7*rdx2SqVol[0]*bcVals[1]+(436425.0*phiUy[0]+4.0567527e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.39120652485041e+7*rdx2SqVolSq[0]*phiUy[1])-1.791344449274544e+7*rdx2SqVolSq[0]*phiLx[1]+3.939936e+7*rdx2SqVolSq[0]*bcVals[1]+(1.630608e+7*phiUy[0]-1.969968e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.813584035008861e+7*rdx2SqVolCu[0]*phiUy[1])-1155172.233549179*rdx2SqVolCu[0]*phiLx[1]+3.9779376e+7*rdx2SqVolCu[0]*bcVals[1]+((-7049385.0*phiUy[0])-1561287.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-8.491586400000001e+7*rdx2SqVolR4[1])+2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.04493464e+8*rdx2SqVolR4[0])*phiPrevC[3]))/(8.491586400000001e+7*rdx2SqVolR4[1]-2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.04493464e+8*rdx2SqVolR4[0]); 

}

void MGpoissonDampedJacobi2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((1.9688856e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+6.447991200000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+(1.840632657003176e+8*rdx2SqVolCu[1]-4.580379954982652e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.214413387187174e+9*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-3.571379299595986e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.601520208069406e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-2.887905122726076e+7*rdx2SqVolCu[0])*rho[1]-1.59900048e+8*rho[0]*rdx2SqVolCu[1]+2.97417736e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.307751256e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+6.7183704e+8*rdx2SqVolCu[0]*rho[0])*omega*volFac+((1.5135219e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+9039240.0*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.9204101e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(8.513594099999999e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.2584276e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.07505661e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+(2.580326154677349e+7*rdx2SqVolR4[1]+2613649.080164168*rdx2SqVol[0]*rdx2SqVolCu[1]-3.383275325638488e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-4.467607425939947e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(9.098581884049787e+7*rdx2SqVol[0]*rdx2SqVolCu[1]-2.41723701273902e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-6.504368232599894e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+(5.243927039999999e+8*rdx2SqVolR4[1]-1.175739312e+9*rdx2SqVol[0]*rdx2SqVolCu[1]-3.78482016e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-8.6379048e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-7448760.0*phiUy[0])-2.54747592e+8*phiPrevC[0])*rdx2SqVolR4[1]+((-1.325858046406458e+7*rdx2SqVol[0]*phiUy[1])-7.395957580471022e+7*rdx2SqVol[0]*phiLx[1]-1.30827312e+8*rdx2SqVol[0]*bcVals[1]+((-4.8756675e+7*phiUy[0])+7.81081488e+8*phiPrevC[0]-7.904150099999999e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-2609403.823634816*rdx2SqVolSq[0]*phiUy[1])+1.629314075353864e+8*rdx2SqVolSq[0]*phiLx[1]+2.03977536e+8*rdx2SqVolSq[0]*bcVals[1]+(2.1840576e+8*phiUy[0]+1.397230704e+9*phiPrevC[0]+1.74784848e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-1.933865037539783e+7*rdx2SqVolCu[0]*phiUy[1])+6.446659452024169e+8*rdx2SqVolCu[0]*phiLx[1]+1.008237552e+9*rdx2SqVolCu[0]*bcVals[1]+(4.49890875e+8*phiUy[0]-1.212108912e+9*phiPrevC[0]+6.899945009999999e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]+3.312597052538735e+8*rdx2SqVolR4[0]*phiLx[1]+5.178539520000001e+8*rdx2SqVolR4[0]*bcVals[1]+(3.54553416e+8*phiLx[0]-6.13480392e+8*phiPrevC[0])*rdx2SqVolR4[0])*omega+2.54747592e+8*phiPrevC[0]*rdx2SqVolR4[1]-7.81081488e+8*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*phiPrevC[0]*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*phiPrevC[0]*rdx2SqVolR4[0])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[1] = (((1.840632657003176e+8*rdx2SqVolCu[1]-1.063183084441229e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-6.255090874156799e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[3]+((-2.6637864e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.723752800000001e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[2]+((-1.59900048e+8*rdx2SqVolCu[1])-3.40569768e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+3.0710148e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.55944656e+8*rdx2SqVolCu[0])*rho[1]+4.831866111218099e+7*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+7.578527340329195e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1]+3.907165754276457e+7*rdx2SqVolCu[0]*rho[0])*omega*volFac+((2.580326154677349e+7*rdx2SqVolR4[1]+2.729876357747648e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.768514841836236e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.03700668718898e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+((-1.307267512076119e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.318918332056357e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4074838.318343807*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[3]+((-2.0477061e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-1.222956e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.5982019e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-1.43100837e+8*rdx2SqVol[0]*rdx2SqVolCu[1])+6.949008e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5507397.0*rdx2SqVolCu[0]*rdx2SqVol[1])*phiLx[2]+((-1.082621267116283e+8*rdx2SqVol[0]*rdx2SqVolCu[1])-2.485380394840879e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-5.023498826926871e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]+((-7448760.0*phiUy[1])-2.54747592e+8*phiPrevC[1])*rdx2SqVolR4[1]+((-2.85606585e+8*rdx2SqVol[0]*phiUy[1])+7.81081488e+8*rdx2SqVol[0]*phiPrevC[1]+1.13565375e+8*rdx2SqVol[0]*phiLx[1]-1.762440955346286e+8*rdx2SqVol[0]*bcVals[1]+(1.793807945138149e+7*phiUy[0]+1.24315031671747e+8*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+(1.7179164e+8*rdx2SqVolSq[0]*phiUy[1]+1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]+2.64231072e+8*rdx2SqVolSq[0]*phiLx[1]-3.358473674432715e+8*rdx2SqVolSq[0]*bcVals[1]+(3530369.879035339*phiUy[0]+2.886623335846444e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.04427225e+8*rdx2SqVolCu[0]*phiUy[1]-1.212108912e+9*rdx2SqVolCu[0]*phiPrevC[1]-1.83058407e+8*rdx2SqVolCu[0]*phiLx[1]+4.004977296097099e+8*rdx2SqVolCu[0]*bcVals[1]+(2.616405639024413e+7*phiUy[0]-2.012954270604647e+8*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1]-6.13480392e+8*rdx2SqVolR4[0]*phiPrevC[1]-9.268408800000001e+7*rdx2SqVolR4[0]*phiLx[1]+2.03852126310076e+8*rdx2SqVolR4[0]*bcVals[1]-1.01926063155038e+8*phiLx[0]*rdx2SqVolR4[0])*omega+2.54747592e+8*phiPrevC[1]*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*phiPrevC[1]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]*phiPrevC[1])/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[2] = -(1.0*(((1.892833619933881e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1043761.529453926*rdx2SqVolSq[0]*rdx2SqVol[1]+2.887905122726076e+7*rdx2SqVolCu[0])*rho[3]+(2979504.0*rdx2SqVolCu[1]+9.3076104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-2.72420232e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-6.7183704e+8*rdx2SqVolCu[0])*rho[2]+((-1.6744728e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-5.4838056e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]-1.565397867170925e+8*rho[0]*rdx2SqVolCu[1]+3.895463326200199e+8*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+1.032818862000306e+9*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2853825.637446513*rdx2SqVol[0]*rdx2SqVolCu[1]+3.600977276616046e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.263458491192659e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(1378128.741702675*rdx2SqVol[0]*rdx2SqVolCu[1]+2.960765571997269e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.267455527520777e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.312597052538735e+8*rdx2SqVolR4[0])*phiLx[3]+(1.02792888e+8*rdx2SqVolR4[1]-2.99762793e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-5.67623952e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.93928705e+8*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0])*phiPrevC[2]+(1472823.0*rdx2SqVol[0]*rdx2SqVolCu[1]+3.1293288e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.35473751e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-3.54553416e+8*rdx2SqVolR4[0])*phiLx[2]+(2.064260923741879e+8*rdx2SqVolR4[1]-3.396327436986036e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.799755648262666e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.163655887686684e+9*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-1.03213046187094e+8*phiUy[0]*rdx2SqVolR4[1]+((-967725.0*rdx2SqVol[0]*phiUy[1])-7.240533299999999e+7*rdx2SqVol[0]*phiLx[1]-1.280780073139847e+8*rdx2SqVol[0]*bcVals[1]+(3.112358382584934e+8*phiUy[0]-7.738046275219913e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.615696e+7*rdx2SqVolSq[0]*phiUy[1])+1.9207188e+8*rdx2SqVolSq[0]*phiLx[1]+3.002634505450664e+8*rdx2SqVolSq[0]*bcVals[1]+(5.441679977753879e+8*phiUy[0]+2.05578101083412e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+(1.5631245e+7*rdx2SqVolCu[0]*phiUy[1]+5.166636929999999e+8*rdx2SqVolCu[0]*phiLx[1]+7.845903330673281e+8*rdx2SqVolCu[0]*bcVals[1]+(5.531752422117667e+8*phiLx[0]-3.636424649020887e+8*phiUy[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-2.54747592e+8*rdx2SqVolR4[1])+7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]-1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]-6.13480392e+8*rdx2SqVolR4[0])*phiPrevC[2]))/(2.54747592e+8*rdx2SqVolR4[1]-7.81081488e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-1.397230704e+9*rdx2SqVolSq[0]*rdx2SqVolSq[1]+1.212108912e+9*rdx2SqVolCu[0]*rdx2SqVol[1]+6.13480392e+8*rdx2SqVolR4[0]); 
  phiC[3] = -(1.0*(((993168.0*rdx2SqVolCu[1]+1.43736648e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-8.4591528e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-5.1981552e+7*rdx2SqVolCu[0])*rho[3]+((-8536308.482054757*rdx2SqVol[0]*rdx2SqVolSq[1])-470715.9838713786*rdx2SqVolSq[0]*rdx2SqVol[1]-1.302388584758819e+7*rdx2SqVolCu[0])*rho[2]+((-5.217992890569751e+7*rdx2SqVolCu[1])+3.014008121001614e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.773250060898033e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*rho[1]+7551544.0*rdx2SqVol[0]*rho[0]*rdx2SqVolSq[1]+2.4730888e+7*rdx2SqVolSq[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((3.4264296e+7*rdx2SqVolR4[1]-8.292744900000001e+7*rdx2SqVol[0]*rdx2SqVolCu[1]+2.5216968e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+2.2741929e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[3]+(8.491586400000001e+7*rdx2SqVolR4[1]-2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.04493464e+8*rdx2SqVolR4[0])*phiPrevC[3]+((-705375.0*rdx2SqVol[0]*rdx2SqVolCu[1])-1.0603404e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5.9861487e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.0894696e+7*rdx2SqVolR4[0])*phiLx[3]+((-1287019.405122937*rdx2SqVol[0]*rdx2SqVolCu[1])-1.623970144356256e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]+5697950.058319833*rdx2SqVolCu[0]*rdx2SqVol[1])*phiUy[2]+((-772143.0538617824*rdx2SqVol[0]*rdx2SqVolCu[1])-1.15968374092867e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+6.553339111300069e+7*rdx2SqVolCu[0]*rdx2SqVol[1]+3.397535438501266e+7*rdx2SqVolR4[0])*phiLx[2]+((-2.4494448e+7*rdx2SqVol[0]*rdx2SqVolCu[1])-3.261216e+7*rdx2SqVolSq[0]*rdx2SqVolSq[1]-2.2558032e+7*rdx2SqVolCu[0]*rdx2SqVol[1])*bcVals[2]-3.440434872903132e+7*phiUy[1]*rdx2SqVolR4[1]+(9.798283298525652e+7*rdx2SqVol[0]*phiUy[1]+3.705960859779652e+7*rdx2SqVol[0]*phiLx[1]-5.7513456e+7*rdx2SqVol[0]*bcVals[1]+(436425.0*phiUy[0]+4.0567527e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolCu[1]+((-3.39120652485041e+7*rdx2SqVolSq[0]*phiUy[1])-1.791344449274544e+7*rdx2SqVolSq[0]*phiLx[1]+3.939936e+7*rdx2SqVolSq[0]*bcVals[1]+(1.630608e+7*phiUy[0]-1.969968e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVolSq[1]+((-2.813584035008861e+7*rdx2SqVolCu[0]*phiUy[1])-1155172.233549179*rdx2SqVolCu[0]*phiLx[1]+3.9779376e+7*rdx2SqVolCu[0]*bcVals[1]+((-7049385.0*phiUy[0])-1561287.0*phiLx[0])*rdx2SqVolCu[0])*rdx2SqVol[1])*omega+((-8.491586400000001e+7*rdx2SqVolR4[1])+2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]+4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]-4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]-2.04493464e+8*rdx2SqVolR4[0])*phiPrevC[3]))/(8.491586400000001e+7*rdx2SqVolR4[1]-2.60360496e+8*rdx2SqVol[0]*rdx2SqVolCu[1]-4.65743568e+8*rdx2SqVolSq[0]*rdx2SqVolSq[1]+4.04036304e+8*rdx2SqVolCu[0]*rdx2SqVol[1]+2.04493464e+8*rdx2SqVolR4[0]); 

}

void MGpoissonDampedJacobi2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])-17043.37994647775*rdx2SqVolSq[0])*rho[1]-92496.0*rho[0]*rdx2SqVolSq[1]-667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(40344.0*phiPrevC[0]-40344.0*phiUy[0])*rdx2SqVolCu[1]+((-50064.92859277839*rdx2SqVol[0]*phiUy[1])-50064.92859277839*rdx2SqVol[0]*phiLx[1]-92496.0*rdx2SqVol[0]*bcVals[1]+((-293355.0*phiUy[0])+345384.0*phiPrevC[0]-52029.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9586.901219893733*rdx2SqVolSq[0]*phiUy[1])-288932.0554646022*rdx2SqVolSq[0]*phiLx[1]-828720.0*rdx2SqVolSq[0]*bcVals[1]+((-52029.0*phiUy[0])+345384.0*phiPrevC[0]-293355.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-39767.88654178142*rdx2SqVolCu[0]*phiLx[1]-115456.0*rdx2SqVolCu[0]*bcVals[1]+(40344.0*phiPrevC[0]-40344.0*phiLx[0])*rdx2SqVolCu[0])*omega-40344.0*phiPrevC[0]*rdx2SqVolCu[1]-345384.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]-345384.0*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]-40344.0*phiPrevC[0]*rdx2SqVolCu[0]))/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-277488.0*rdx2SqVolSq[1])-426384.0*rdx2SqVol[0]*rdx2SqVol[1]-53136.0*rdx2SqVolSq[0])*rho[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(121032.0*phiPrevC[1]-121032.0*phiUy[1])*rdx2SqVolCu[1]+((-221031.0*rdx2SqVol[0]*phiUy[1])+1036152.0*rdx2SqVol[0]*phiPrevC[1]+167649.0*rdx2SqVol[0]*phiLx[1]-373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-29889.0*rdx2SqVolSq[0]*phiUy[1])+1036152.0*rdx2SqVolSq[0]*phiPrevC[1]+11367.0*rdx2SqVolSq[0]*phiLx[1]-1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiUy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]*phiPrevC[1]-2952.0*rdx2SqVolCu[0]*phiLx[1]-136347.039571822*rdx2SqVolCu[0]*bcVals[1])*omega-121032.0*phiPrevC[1]*rdx2SqVolCu[1]-1036152.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0]*phiPrevC[1]))/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]-95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-64764.84379661545*rho[0]*rdx2SqVolSq[1]-454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+((-35055.0*rdx2SqVol[0]*phiUy[1])-35055.0*rdx2SqVol[0]*phiLx[1]-64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiUy[0]-36430.22463559618*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(35055.0*rdx2SqVolSq[0]*phiUy[1]-188385.0*rdx2SqVolSq[0]*phiLx[1]-583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiUy[0]-190246.7286525579*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[2])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+((-64764.84379661545*rdx2SqVolSq[1])-34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]-121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiUy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[2]+(190246.7286525579*rdx2SqVol[0]*phiUy[1]+39128.75979378851*rdx2SqVol[0]*phiLx[1]-87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiUy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(36430.22463559618*rdx2SqVolSq[0]*phiUy[1]-44657.46597154836*rdx2SqVolSq[0]*phiLx[1]-168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiUy[0]-44403.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((980220.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-378861.8654443858*rdx2SqVolSq[1])-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2012649.966598266*rdx2SqVol[0]*rdx2SqVol[1])-378861.8654443858*rdx2SqVolSq[0])*rho[1]+924336.0*rho[0]*rdx2SqVolSq[1]+5185588.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+924336.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1381887.0*rdx2SqVol[0]*rdx2SqVolSq[1])-41013.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-41013.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1381887.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1862784.0*rdx2SqVolCu[1]+1.1134104e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+4159512.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(549960.7724192695*rdx2SqVolCu[1]+2951378.203030407*rdx2SqVol[0]*rdx2SqVolSq[1]+100062.3072040616*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-71036.59977082233*rdx2SqVol[0]*rdx2SqVolSq[1])-1544603.07302135*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(624456.0*phiLy[0]-1555848.0*phiPrevC[0])*rdx2SqVolCu[1]+((-1544603.07302135*rdx2SqVol[0]*phiLy[1])+100062.3072040616*rdx2SqVol[0]*phiLx[1]+4159512.0*rdx2SqVol[0]*bcVals[1]+((-1.1189052e+7*phiPrevC[0])+3368931.0*phiLy[0]+173313.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-71036.59977082233*rdx2SqVolSq[0]*phiLy[1])+2951378.203030407*rdx2SqVolSq[0]*phiLx[1]+1.1134104e+7*rdx2SqVolSq[0]*bcVals[1]+((-1.1189052e+7*phiPrevC[0])+173313.0*phiLy[0]+3368931.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+549960.7724192695*rdx2SqVolCu[0]*phiLx[1]+1862784.0*rdx2SqVolCu[0]*bcVals[1]+(624456.0*phiLx[0]-1555848.0*phiPrevC[0])*rdx2SqVolCu[0])*omega+1555848.0*phiPrevC[0]*rdx2SqVolCu[1]+1.1189052e+7*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*phiPrevC[0]*rdx2SqVolCu[0])/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((378861.8654443858*rdx2SqVolSq[1]+333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-924336.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-275184.0*rdx2SqVolSq[0])*rho[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+207762.9584695019*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-549960.7724192695*rdx2SqVolCu[1])-583616.2516611406*rdx2SqVol[0]*rdx2SqVolSq[1]-29789.54183937711*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-449898.4652152082*rdx2SqVol[0]*rdx2SqVolSq[1])-453764.4026177019*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1708037.655172742*rdx2SqVol[0]*rdx2SqVolSq[1]+934933.3131127583*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(757809.0*rdx2SqVol[0]*rdx2SqVolSq[1]+22491.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-451143.0*rdx2SqVol[0]*rdx2SqVolSq[1])-497457.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(1555848.0*phiPrevC[1]-624456.0*phiLy[1])*rdx2SqVolCu[1]+(1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]-722367.0*rdx2SqVol[0]*phiLy[1]+1097649.0*rdx2SqVol[0]*phiLx[1]-5603489.203427449*rdx2SqVol[0]*bcVals[1]+(847040.3948826761*phiLy[0]+1100685.379244677*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]-51597.0*rdx2SqVolSq[0]*phiLy[1]+2182239.0*rdx2SqVolSq[0]*phiLx[1]-5563665.891259826*rdx2SqVolSq[0]*bcVals[1]+(38955.5547130316*phiLy[0]+2275410.734360502*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]*phiPrevC[1]+349272.0*rdx2SqVolCu[0]*phiLx[1]-733281.0298923596*rdx2SqVolCu[0]*bcVals[1]+366640.5149461798*phiLx[0]*rdx2SqVolCu[0])*omega-1555848.0*phiPrevC[1]*rdx2SqVolCu[1]-1.1189052e+7*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0]*phiPrevC[1]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((333114.9395148767*rdx2SqVol[0]*rdx2SqVol[1]+378861.8654443858*rdx2SqVolSq[0])*rho[3]+((-275184.0*rdx2SqVolSq[1])-1737060.0*rdx2SqVol[0]*rdx2SqVol[1]-924336.0*rdx2SqVolSq[0])*rho[2]-537540.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+207762.9584695019*rho[0]*rdx2SqVolSq[1]+1103711.272005501*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-453764.4026177019*rdx2SqVol[0]*rdx2SqVolSq[1])-449898.4652152082*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-29789.54183937711*rdx2SqVol[0]*rdx2SqVolSq[1])-583616.2516611406*rdx2SqVolSq[0]*rdx2SqVol[1]-549960.7724192695*rdx2SqVolCu[0])*phiLx[3]+((-733281.0298923596*rdx2SqVolCu[1])-5563665.891259826*rdx2SqVol[0]*rdx2SqVolSq[1]-5603489.203427449*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0])*phiPrevC[2]+(349272.0*rdx2SqVolCu[1]+2182239.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1097649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-51597.0*rdx2SqVol[0]*rdx2SqVolSq[1])-722367.0*rdx2SqVolSq[0]*rdx2SqVol[1]-624456.0*rdx2SqVolCu[0])*phiLx[2]+366640.5149461798*phiLy[0]*rdx2SqVolCu[1]+((-497457.0*rdx2SqVol[0]*phiLy[1])+22491.0*rdx2SqVol[0]*phiLx[1]+934933.3131127583*rdx2SqVol[0]*bcVals[1]+(2275410.734360502*phiLy[0]+38955.5547130316*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-451143.0*rdx2SqVolSq[0]*phiLy[1])+757809.0*rdx2SqVolSq[0]*phiLx[1]+1708037.655172742*rdx2SqVolSq[0]*bcVals[1]+(1100685.379244677*phiLy[0]+847040.3948826761*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1555848.0*rdx2SqVolCu[1])-1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-1555848.0*rdx2SqVolCu[0])*phiPrevC[2]))/(1555848.0*rdx2SqVolCu[1]+1.1189052e+7*rdx2SqVol[0]*rdx2SqVolSq[1]+1.1189052e+7*rdx2SqVolSq[0]*rdx2SqVol[1]+1555848.0*rdx2SqVolCu[0]); 
  phiC[3] = (((91728.0*rdx2SqVolSq[1]+388764.0*rdx2SqVol[0]*rdx2SqVol[1]+91728.0*rdx2SqVolSq[0])*rho[3]+((-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])-69254.31948983399*rdx2SqVolSq[0])*rho[2]+((-69254.31948983399*rdx2SqVolSq[1])-60891.97819089144*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+98260.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-518616.0*rdx2SqVolCu[1])-3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]-3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]-518616.0*rdx2SqVolCu[0])*phiPrevC[3]+((-116424.0*rdx2SqVolCu[1])-468249.0*rdx2SqVol[0]*rdx2SqVolSq[1]-108927.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-108927.0*rdx2SqVol[0]*rdx2SqVolSq[1])-468249.0*rdx2SqVolSq[0]*rdx2SqVol[1]-116424.0*rdx2SqVolCu[0])*phiLx[3]+(73032.0*rdx2SqVol[0]*rdx2SqVolSq[1]-419832.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(82946.18112366594*rdx2SqVol[0]*rdx2SqVolSq[1]+82239.50439417786*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-109228.3200777161*rdx2SqVol[0]*rdx2SqVolSq[1])-474351.5585164657*rdx2SqVolSq[0]*rdx2SqVol[1]-122213.5049820599*rdx2SqVolCu[0])*phiLx[2]-122213.5049820599*phiLy[1]*rdx2SqVolCu[1]+((-474351.5585164657*rdx2SqVol[0]*phiLy[1])+82239.50439417786*rdx2SqVol[0]*phiLx[1]-419832.0*rdx2SqVol[0]*bcVals[1]+(90933.0*phiLy[0]+82467.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-109228.3200777161*rdx2SqVolSq[0]*phiLy[1])+82946.18112366594*rdx2SqVolSq[0]*phiLx[1]+73032.0*rdx2SqVolSq[0]*bcVals[1]+(82467.0*phiLy[0]+90933.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0])*phiPrevC[3])/(518616.0*rdx2SqVolCu[1]+3729684.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3729684.0*rdx2SqVolSq[0]*rdx2SqVol[1]+518616.0*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1463496.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-6483052.316323845*rdx2SqVolSq[1])-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])-6483052.316323845*rdx2SqVolSq[0])*rho[1]+1.5082056e+8*rho[0]*rdx2SqVolSq[1]+6.76239224e+8*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+1.5082056e+8*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1.3788513e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-2998647.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-2998647.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1.3788513e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1.16252928e+8*rdx2SqVolCu[1]+5.22905808e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.2339864e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(7.436442362842056e+7*rdx2SqVolCu[1]+3.323615026097367e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+6.975998306284967e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-3204690.633637355*rdx2SqVol[0]*rdx2SqVolSq[1])-1.476291855196238e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(7.959362400000001e+7*phiLy[0]-1.37720088e+8*phiPrevC[0])*rdx2SqVolCu[1]+((-1.476291855196238e+7*rdx2SqVol[0]*phiLy[1])+6.975998306284967e+7*rdx2SqVol[0]*phiLx[1]+1.2339864e+8*rdx2SqVol[0]*bcVals[1]+((-7.53412248e+8*phiPrevC[0])+3.55706679e+8*phiLy[0]+7.4553345e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-3204690.633637355*rdx2SqVolSq[0]*phiLy[1])+3.323615026097367e+8*rdx2SqVolSq[0]*phiLx[1]+5.22905808e+8*rdx2SqVolSq[0]*bcVals[1]+((-7.53412248e+8*phiPrevC[0])+7.4553345e+7*phiLy[0]+3.55706679e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+7.436442362842056e+7*rdx2SqVolCu[0]*phiLx[1]+1.16252928e+8*rdx2SqVolCu[0]*bcVals[1]+(7.959362400000001e+7*phiLx[0]-1.37720088e+8*phiPrevC[0])*rdx2SqVolCu[0])*omega+1.37720088e+8*phiPrevC[0]*rdx2SqVolCu[1]+7.53412248e+8*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*phiPrevC[0]*rdx2SqVolCu[0])/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((6483052.316323845*rdx2SqVolSq[1]+1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-1.5082056e+8*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-3.5007984e+7*rdx2SqVolSq[0])*rho[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-8771188.427967556*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-7.436442362842056e+7*rdx2SqVolCu[1])-8.604493260170916e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.619246322188773e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-4604440.565570912*rdx2SqVol[0]*rdx2SqVolSq[1])-92486.31697175531*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+((-2.832900731640274e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-7176426.895609816*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-1.8655047e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-4056993.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-5040279.0*rdx2SqVol[0]*rdx2SqVolSq[1])-125001.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(1.37720088e+8*phiPrevC[1]-7.959362400000001e+7*phiLy[1])*rdx2SqVolCu[1]+(7.53412248e+8*rdx2SqVol[0]*phiPrevC[1]-9.1983429e+7*rdx2SqVol[0]*phiLy[1]+1.07116875e+8*rdx2SqVol[0]*phiLx[1]-1.662365553838119e+8*rdx2SqVol[0]*bcVals[1]+(1.172561417439388e+8*phiLx[0]-1.997336039383145e+7*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(7.53412248e+8*rdx2SqVolSq[0]*phiPrevC[1]-1.7305083e+7*rdx2SqVolSq[0]*phiLy[1]+1.13325453e+8*rdx2SqVolSq[0]*phiLx[1]-2.331518580374791e+8*rdx2SqVolSq[0]*bcVals[1]+(1.244999003826421e+8*phiLx[0]-4335757.916097598*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]*phiPrevC[1]+2.0806632e+7*rdx2SqVolCu[0]*phiLx[1]-4.576272223287419e+7*rdx2SqVolCu[0]*bcVals[1]+2.28813611164371e+7*phiLx[0]*rdx2SqVolCu[0])*omega-1.37720088e+8*phiPrevC[1]*rdx2SqVolCu[1]-7.53412248e+8*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-7.53412248e+8*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]-1.37720088e+8*rdx2SqVolCu[0]*phiPrevC[1]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1]+6483052.316323845*rdx2SqVolSq[0])*rho[3]+((-3.5007984e+7*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-1.5082056e+8*rdx2SqVolSq[0])*rho[2]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-8771188.427967556*rho[0]*rdx2SqVolSq[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-92486.31697175531*rdx2SqVol[0]*rdx2SqVolSq[1])-4604440.565570912*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1.619246322188773e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.604493260170916e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.436442362842056e+7*rdx2SqVolCu[0])*phiLx[3]+((-4.576272223287419e+7*rdx2SqVolCu[1])-2.331518580374791e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.662365553838119e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0])*phiPrevC[2]+(2.0806632e+7*rdx2SqVolCu[1]+1.13325453e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.07116875e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-1.7305083e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-9.1983429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.959362400000001e+7*rdx2SqVolCu[0])*phiLx[2]+2.28813611164371e+7*phiLy[0]*rdx2SqVolCu[1]+((-125001.0*rdx2SqVol[0]*phiLy[1])-4056993.0*rdx2SqVol[0]*phiLx[1]-7176426.895609816*rdx2SqVol[0]*bcVals[1]+(1.244999003826421e+8*phiLy[0]-4335757.916097598*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5040279.0*rdx2SqVolSq[0]*phiLy[1])-1.8655047e+7*rdx2SqVolSq[0]*phiLx[1]-2.832900731640274e+7*rdx2SqVolSq[0]*bcVals[1]+(1.172561417439388e+8*phiLy[0]-1.997336039383145e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1.37720088e+8*rdx2SqVolCu[1])-7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-1.37720088e+8*rdx2SqVolCu[0])*phiPrevC[2]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[3] = (((1.1669328e+7*rdx2SqVolSq[1]+5.2828968e+7*rdx2SqVol[0]*rdx2SqVol[1]+1.1669328e+7*rdx2SqVolSq[0])*rho[3]+(640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1]+2923729.475989186*rdx2SqVolSq[0])*rho[2]+(2923729.475989186*rdx2SqVolSq[1]+640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+892952.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-4.5906696e+7*rdx2SqVolCu[1])-2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-4.5906696e+7*rdx2SqVolCu[0])*phiPrevC[3]+((-6935544.0*rdx2SqVolCu[1])-3.7224429e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-8287875.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-8287875.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3.7224429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-6935544.0*rdx2SqVolCu[0])*phiLx[3]+(1436304.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3222576.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-41709.51549706613*rdx2SqVol[0]*rdx2SqVolSq[1])-2076512.411924138*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-9072373.0108449*rdx2SqVol[0]*rdx2SqVolSq[1])-4.075563181606774e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7627120.3721457*rdx2SqVolCu[0])*phiLx[2]-7627120.3721457*phiLy[1]*rdx2SqVolCu[1]+((-4.075563181606774e+7*rdx2SqVol[0]*phiLy[1])-2076512.411924138*rdx2SqVol[0]*phiLx[1]+3222576.0*rdx2SqVol[0]*bcVals[1]+((-56373.0*phiLy[0])-2273067.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9072373.0108449*rdx2SqVolSq[0]*phiLy[1])-41709.51549706613*rdx2SqVolSq[0]*phiLx[1]+1436304.0*rdx2SqVolSq[0]*bcVals[1]+((-2273067.0*phiLy[0])-56373.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(4.5906696e+7*rdx2SqVolCu[1]+2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+4.5906696e+7*rdx2SqVolCu[0])*phiPrevC[3])/(4.5906696e+7*rdx2SqVolCu[1]+2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+4.5906696e+7*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1463496.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+((-6483052.316323845*rdx2SqVolSq[1])-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+((-2.756345471586376e+7*rdx2SqVol[0]*rdx2SqVol[1])-6483052.316323845*rdx2SqVolSq[0])*rho[1]+1.5082056e+8*rho[0]*rdx2SqVolSq[1]+6.76239224e+8*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+1.5082056e+8*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-1.3788513e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-2998647.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-2998647.0*rdx2SqVol[0]*rdx2SqVolSq[1])-1.3788513e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(1.16252928e+8*rdx2SqVolCu[1]+5.22905808e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.2339864e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(7.436442362842056e+7*rdx2SqVolCu[1]+3.323615026097367e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+6.975998306284967e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-3204690.633637355*rdx2SqVol[0]*rdx2SqVolSq[1])-1.476291855196238e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(7.959362400000001e+7*phiLy[0]-1.37720088e+8*phiPrevC[0])*rdx2SqVolCu[1]+((-1.476291855196238e+7*rdx2SqVol[0]*phiLy[1])+6.975998306284967e+7*rdx2SqVol[0]*phiLx[1]+1.2339864e+8*rdx2SqVol[0]*bcVals[1]+((-7.53412248e+8*phiPrevC[0])+3.55706679e+8*phiLy[0]+7.4553345e+7*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-3204690.633637355*rdx2SqVolSq[0]*phiLy[1])+3.323615026097367e+8*rdx2SqVolSq[0]*phiLx[1]+5.22905808e+8*rdx2SqVolSq[0]*bcVals[1]+((-7.53412248e+8*phiPrevC[0])+7.4553345e+7*phiLy[0]+3.55706679e+8*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+7.436442362842056e+7*rdx2SqVolCu[0]*phiLx[1]+1.16252928e+8*rdx2SqVolCu[0]*bcVals[1]+(7.959362400000001e+7*phiLx[0]-1.37720088e+8*phiPrevC[0])*rdx2SqVolCu[0])*omega+1.37720088e+8*phiPrevC[0]*rdx2SqVolCu[1]+7.53412248e+8*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*phiPrevC[0]*rdx2SqVolCu[0])/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[1] = -(1.0*(((6483052.316323845*rdx2SqVolSq[1]+1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+((-1.5082056e+8*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-3.5007984e+7*rdx2SqVolSq[0])*rho[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]-8771188.427967556*rdx2SqVolSq[0]*rho[0])*omega*volFac+(((-7.436442362842056e+7*rdx2SqVolCu[1])-8.604493260170916e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-1.619246322188773e+7*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-4604440.565570912*rdx2SqVol[0]*rdx2SqVolSq[1])-92486.31697175531*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+((-2.832900731640274e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-7176426.895609816*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-1.8655047e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-4056993.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-5040279.0*rdx2SqVol[0]*rdx2SqVolSq[1])-125001.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(1.37720088e+8*phiPrevC[1]-7.959362400000001e+7*phiLy[1])*rdx2SqVolCu[1]+(7.53412248e+8*rdx2SqVol[0]*phiPrevC[1]-9.1983429e+7*rdx2SqVol[0]*phiLy[1]+1.07116875e+8*rdx2SqVol[0]*phiLx[1]-1.662365553838119e+8*rdx2SqVol[0]*bcVals[1]+(1.172561417439388e+8*phiLx[0]-1.997336039383145e+7*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(7.53412248e+8*rdx2SqVolSq[0]*phiPrevC[1]-1.7305083e+7*rdx2SqVolSq[0]*phiLy[1]+1.13325453e+8*rdx2SqVolSq[0]*phiLx[1]-2.331518580374791e+8*rdx2SqVolSq[0]*bcVals[1]+(1.244999003826421e+8*phiLx[0]-4335757.916097598*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]*phiPrevC[1]+2.0806632e+7*rdx2SqVolCu[0]*phiLx[1]-4.576272223287419e+7*rdx2SqVolCu[0]*bcVals[1]+2.28813611164371e+7*phiLx[0]*rdx2SqVolCu[0])*omega-1.37720088e+8*phiPrevC[1]*rdx2SqVolCu[1]-7.53412248e+8*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]-7.53412248e+8*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]-1.37720088e+8*rdx2SqVolCu[0]*phiPrevC[1]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[2] = -(1.0*(((1419713.549541597*rdx2SqVol[0]*rdx2SqVol[1]+6483052.316323845*rdx2SqVolSq[0])*rho[3]+((-3.5007984e+7*rdx2SqVolSq[1])-1.8384852e+8*rdx2SqVol[0]*rdx2SqVol[1]-1.5082056e+8*rdx2SqVolSq[0])*rho[2]+1980024.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]-8771188.427967556*rho[0]*rdx2SqVolSq[1]-3.729173285087449e+7*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-92486.31697175531*rdx2SqVol[0]*rdx2SqVolSq[1])-4604440.565570912*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-1.619246322188773e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-8.604493260170916e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.436442362842056e+7*rdx2SqVolCu[0])*phiLx[3]+((-4.576272223287419e+7*rdx2SqVolCu[1])-2.331518580374791e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-1.662365553838119e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0])*phiPrevC[2]+(2.0806632e+7*rdx2SqVolCu[1]+1.13325453e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+1.07116875e+8*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-1.7305083e+7*rdx2SqVol[0]*rdx2SqVolSq[1])-9.1983429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7.959362400000001e+7*rdx2SqVolCu[0])*phiLx[2]+2.28813611164371e+7*phiLy[0]*rdx2SqVolCu[1]+((-125001.0*rdx2SqVol[0]*phiLy[1])-4056993.0*rdx2SqVol[0]*phiLx[1]-7176426.895609816*rdx2SqVol[0]*bcVals[1]+(1.244999003826421e+8*phiLy[0]-4335757.916097598*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-5040279.0*rdx2SqVolSq[0]*phiLy[1])-1.8655047e+7*rdx2SqVolSq[0]*phiLx[1]-2.832900731640274e+7*rdx2SqVolSq[0]*bcVals[1]+(1.172561417439388e+8*phiLy[0]-1.997336039383145e+7*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+((-1.37720088e+8*rdx2SqVolCu[1])-7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-1.37720088e+8*rdx2SqVolCu[0])*phiPrevC[2]))/(1.37720088e+8*rdx2SqVolCu[1]+7.53412248e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+7.53412248e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+1.37720088e+8*rdx2SqVolCu[0]); 
  phiC[3] = (((1.1669328e+7*rdx2SqVolSq[1]+5.2828968e+7*rdx2SqVol[0]*rdx2SqVol[1]+1.1669328e+7*rdx2SqVolSq[0])*rho[3]+(640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1]+2923729.475989186*rdx2SqVolSq[0])*rho[2]+(2923729.475989186*rdx2SqVolSq[1]+640262.9733226809*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+892952.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-4.5906696e+7*rdx2SqVolCu[1])-2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]-2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]-4.5906696e+7*rdx2SqVolCu[0])*phiPrevC[3]+((-6935544.0*rdx2SqVolCu[1])-3.7224429e+7*rdx2SqVol[0]*rdx2SqVolSq[1]-8287875.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-8287875.0*rdx2SqVol[0]*rdx2SqVolSq[1])-3.7224429e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-6935544.0*rdx2SqVolCu[0])*phiLx[3]+(1436304.0*rdx2SqVol[0]*rdx2SqVolSq[1]+3222576.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-41709.51549706613*rdx2SqVol[0]*rdx2SqVolSq[1])-2076512.411924138*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-9072373.0108449*rdx2SqVol[0]*rdx2SqVolSq[1])-4.075563181606774e+7*rdx2SqVolSq[0]*rdx2SqVol[1]-7627120.3721457*rdx2SqVolCu[0])*phiLx[2]-7627120.3721457*phiLy[1]*rdx2SqVolCu[1]+((-4.075563181606774e+7*rdx2SqVol[0]*phiLy[1])-2076512.411924138*rdx2SqVol[0]*phiLx[1]+3222576.0*rdx2SqVol[0]*bcVals[1]+((-56373.0*phiLy[0])-2273067.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-9072373.0108449*rdx2SqVolSq[0]*phiLy[1])-41709.51549706613*rdx2SqVolSq[0]*phiLx[1]+1436304.0*rdx2SqVolSq[0]*bcVals[1]+((-2273067.0*phiLy[0])-56373.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(4.5906696e+7*rdx2SqVolCu[1]+2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+4.5906696e+7*rdx2SqVolCu[0])*phiPrevC[3])/(4.5906696e+7*rdx2SqVolCu[1]+2.511374160000001e+8*rdx2SqVol[0]*rdx2SqVolSq[1]+2.511374160000001e+8*rdx2SqVolSq[0]*rdx2SqVol[1]+4.5906696e+7*rdx2SqVolCu[0]); 

}

void MGpoissonDampedJacobi2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  double rdx2SqVolCu[2]; 
  double rdx2SqVolR4[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVolCu[0] = rdx2SqVol[0]*rdx2SqVolSq[0]; 
  rdx2SqVolR4[0] = rdx2SqVol[0]*rdx2SqVolCu[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 
  rdx2SqVolCu[1] = rdx2SqVol[1]*rdx2SqVolSq[1]; 
  rdx2SqVolR4[1] = rdx2SqVol[1]*rdx2SqVolCu[1]; 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((25200.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[3]+(17043.37994647775*rdx2SqVolSq[1]+119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1])*rho[2]+(119719.3518191607*rdx2SqVol[0]*rdx2SqVol[1]+17043.37994647775*rdx2SqVolSq[0])*rho[1]+92496.0*rho[0]*rdx2SqVolSq[1]+667440.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+92496.0*rdx2SqVolSq[0]*rho[0])*omega*volFac+((49575.0*rdx2SqVol[0]*rdx2SqVolSq[1]+9225.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(9225.0*rdx2SqVol[0]*rdx2SqVolSq[1]+49575.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[3]+(115456.0*rdx2SqVolCu[1]+828720.0*rdx2SqVol[0]*rdx2SqVolSq[1]+92496.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(39767.88654178142*rdx2SqVolCu[1]+288932.0554646022*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(9586.901219893733*rdx2SqVol[0]*rdx2SqVolSq[1]+50064.92859277839*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+(40344.0*phiLy[0]-40344.0*phiPrevC[0])*rdx2SqVolCu[1]+(50064.92859277839*rdx2SqVol[0]*phiLy[1]+50064.92859277839*rdx2SqVol[0]*phiLx[1]+92496.0*rdx2SqVol[0]*bcVals[1]+((-345384.0*phiPrevC[0])+293355.0*phiLy[0]+52029.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+(9586.901219893733*rdx2SqVolSq[0]*phiLy[1]+288932.0554646022*rdx2SqVolSq[0]*phiLx[1]+828720.0*rdx2SqVolSq[0]*bcVals[1]+((-345384.0*phiPrevC[0])+52029.0*phiLy[0]+293355.0*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]+39767.88654178142*rdx2SqVolCu[0]*phiLx[1]+115456.0*rdx2SqVolCu[0]*bcVals[1]+(40344.0*phiLx[0]-40344.0*phiPrevC[0])*rdx2SqVolCu[0])*omega+40344.0*phiPrevC[0]*rdx2SqVolCu[1]+345384.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*phiPrevC[0]*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*phiPrevC[0]*rdx2SqVolCu[0])/(40344.0*rdx2SqVolCu[1]+345384.0*rdx2SqVol[0]*rdx2SqVolSq[1]+345384.0*rdx2SqVolSq[0]*rdx2SqVol[1]+40344.0*rdx2SqVolCu[0]); 
  phiC[1] = (((51130.13983943324*rdx2SqVolSq[1]+27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1])*rho[3]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[2]+(277488.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0]*rho[0])*omega*volFac+((119303.6596253443*rdx2SqVolCu[1]+214211.3836260809*rdx2SqVol[0]*rdx2SqVolSq[1]+28760.70365968119*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(35255.89418806449*rdx2SqVolSq[0]*rdx2SqVol[1]-30891.12615299092*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[3]+(583936.681060541*rdx2SqVol[0]*rdx2SqVolSq[1]+64764.84379661545*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(188385.0*rdx2SqVol[0]*rdx2SqVolSq[1]+35055.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(35055.0*rdx2SqVolSq[0]*rdx2SqVol[1]-35055.0*rdx2SqVol[0]*rdx2SqVolSq[1])*phiLx[2]+(121032.0*phiLy[1]-121032.0*phiPrevC[1])*rdx2SqVolCu[1]+((-1036152.0*rdx2SqVol[0]*phiPrevC[1])+221031.0*rdx2SqVol[0]*phiLy[1]-167649.0*rdx2SqVol[0]*phiLx[1]+373818.1334927453*rdx2SqVol[0]*bcVals[1]+(190246.7286525579*phiLy[0]-190246.7286525579*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-1036152.0*rdx2SqVolSq[0]*phiPrevC[1])+29889.0*rdx2SqVolSq[0]*phiLy[1]-11367.0*rdx2SqVolSq[0]*phiLx[1]+1029337.010328493*rdx2SqVolSq[0]*bcVals[1]+(36430.22463559618*phiLy[0]-36430.22463559618*phiLx[0])*rdx2SqVolSq[0])*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0]*phiPrevC[1]+2952.0*rdx2SqVolCu[0]*phiLx[1]+136347.039571822*rdx2SqVolCu[0]*bcVals[1])*omega+121032.0*phiPrevC[1]*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*phiPrevC[1]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]*phiPrevC[1])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[2] = (((27435.68479189101*rdx2SqVol[0]*rdx2SqVol[1]+51130.13983943324*rdx2SqVolSq[0])*rho[3]+(53136.0*rdx2SqVolSq[1]+426384.0*rdx2SqVol[0]*rdx2SqVol[1]+277488.0*rdx2SqVolSq[0])*rho[2]+95760.0*rdx2SqVol[0]*rdx2SqVol[1]*rho[1]+64764.84379661545*rho[0]*rdx2SqVolSq[1]+454933.5369128108*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+((35255.89418806449*rdx2SqVol[0]*rdx2SqVolSq[1]-30891.12615299092*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+(28760.70365968119*rdx2SqVol[0]*rdx2SqVolSq[1]+214211.3836260809*rdx2SqVolSq[0]*rdx2SqVol[1]+119303.6596253443*rdx2SqVolCu[0])*phiLx[3]+(136347.039571822*rdx2SqVolCu[1]+1029337.010328493*rdx2SqVol[0]*rdx2SqVolSq[1]+373818.1334927453*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[2]+(2952.0*rdx2SqVolCu[1]-11367.0*rdx2SqVol[0]*rdx2SqVolSq[1]-167649.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+(29889.0*rdx2SqVol[0]*rdx2SqVolSq[1]+221031.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiLx[2]+(35055.0*rdx2SqVol[0]*phiLy[1]+35055.0*rdx2SqVol[0]*phiLx[1]+64764.84379661545*rdx2SqVol[0]*bcVals[1]+(36430.22463559618*phiLx[0]-36430.22463559618*phiLy[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-35055.0*rdx2SqVolSq[0]*phiLy[1])+188385.0*rdx2SqVolSq[0]*phiLx[1]+583936.681060541*rdx2SqVolSq[0]*bcVals[1]+(190246.7286525579*phiLx[0]-190246.7286525579*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[2])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 
  phiC[3] = (((53136.0*rdx2SqVolSq[1]+306000.0*rdx2SqVol[0]*rdx2SqVol[1]+53136.0*rdx2SqVolSq[0])*rho[3]+(34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1]+64764.84379661545*rdx2SqVolSq[0])*rho[2]+(64764.84379661545*rdx2SqVolSq[1]+34751.86740306195*rdx2SqVol[0]*rdx2SqVol[1])*rho[1]+121296.0*rdx2SqVol[0]*rho[0]*rdx2SqVol[1])*omega*volFac+(((-121032.0*rdx2SqVolCu[1])-1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]-1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]-121032.0*rdx2SqVolCu[0])*phiPrevC[3]+(2952.0*rdx2SqVolCu[1]-166065.0*rdx2SqVol[0]*rdx2SqVolSq[1]-32103.0*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[3]+((-32103.0*rdx2SqVol[0]*rdx2SqVolSq[1])-166065.0*rdx2SqVolSq[0]*rdx2SqVol[1]+2952.0*rdx2SqVolCu[0])*phiLx[3]+(168112.0*rdx2SqVol[0]*rdx2SqVolSq[1]+87248.0*rdx2SqVolSq[0]*rdx2SqVol[1])*bcVals[3]+(44657.46597154836*rdx2SqVol[0]*rdx2SqVolSq[1]-39128.75979378851*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLy[2]+((-36430.22463559618*rdx2SqVol[0]*rdx2SqVolSq[1])-190246.7286525579*rdx2SqVolSq[0]*rdx2SqVol[1])*phiLx[2]+((-190246.7286525579*rdx2SqVol[0]*phiLy[1])-39128.75979378851*rdx2SqVol[0]*phiLx[1]+87248.0*rdx2SqVol[0]*bcVals[1]+(44403.0*phiLy[0]-44403.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVolSq[1]+((-36430.22463559618*rdx2SqVolSq[0]*phiLy[1])+44657.46597154836*rdx2SqVolSq[0]*phiLx[1]+168112.0*rdx2SqVolSq[0]*bcVals[1]+(44403.0*phiLx[0]-44403.0*phiLy[0])*rdx2SqVolSq[0])*rdx2SqVol[1])*omega+(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0])*phiPrevC[3])/(121032.0*rdx2SqVolCu[1]+1036152.0*rdx2SqVol[0]*rdx2SqVolSq[1]+1036152.0*rdx2SqVolSq[0]*rdx2SqVol[1]+121032.0*rdx2SqVolCu[0]); 

}

void MGpoissonResidue2xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1]-18.0*phiC[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]+8.660254037844386*rdx2SqVol[0]*phiLx[3]-7.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]-7.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-46.0*rdx2SqVol[1])-18.0*rdx2SqVol[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdx2SqVol[1]*phiUy[3]-7.0*rdx2SqVol[0]*phiUx[3]-7.0*rdx2SqVol[1]*phiLy[3]-7.0*rdx2SqVol[0]*phiLx[3]+((-46.0*rdx2SqVol[1])-46.0*rdx2SqVol[0])*phiC[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdx2SqVol[1]-1.732050807568877*rdx2SqVol[0]*phiUx[1]+53.6935750346352*rdx2SqVol[0]*phiC[1]+(3.0*phiUx[0]-39.0*phiC[0]+72.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1]-18.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiUx[1]-131.0*rdx2SqVol[0]*phiC[1]+(19.05255888325765*phiUx[0]+29.44486372867091*phiC[0]-96.99484522385713*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-1.732050807568877*rdx2SqVol[0]*phiUx[3]+53.6935750346352*rdx2SqVol[0]*phiC[3]-7.0*rdx2SqVol[1]*phiUy[2]+3.0*rdx2SqVol[0]*phiUx[2]-7.0*rdx2SqVol[1]*phiLy[2]+((-46.0*rdx2SqVol[1])-39.0*rdx2SqVol[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdx2SqVol[1]*phiUy[3]-19.0*rdx2SqVol[0]*phiUx[3]-7.0*rdx2SqVol[1]*phiLy[3]+((-46.0*rdx2SqVol[1])-131.0*rdx2SqVol[0])*phiC[3]+19.05255888325765*rdx2SqVol[0]*phiUx[2]+29.44486372867091*rdx2SqVol[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]-8.660254037844386*rdx2SqVol[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac-77.94228634059945*rdx2SqVol[1]*phiUy[3]+77.94228634059945*rdx2SqVol[1]*phiLy[3]+(81.0*phiUy[1]+81.0*phiLy[1]-162.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiUx[1]-423.0*rdx2SqVol[0]*phiC[1]+(98.726896031426*phiUx[0]-98.726896031426*phiC[0]+193.9896904477143*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]-8.660254037844386*rdx2SqVol[0]*phiC[3]-7.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]-7.0*rdx2SqVol[1]*phiLy[2]+((-46.0*rdx2SqVol[1])-9.0*rdx2SqVol[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-21.0*rdx2SqVol[1]*phiUy[3]-29.0*rdx2SqVol[0]*phiUx[3]-21.0*rdx2SqVol[1]*phiLy[3]+((-138.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]+32.90896534380867*rdx2SqVol[0]*phiUx[2]-32.90896534380867*rdx2SqVol[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdx2SqVol[1]+1.732050807568877*rdx2SqVol[0]*phiLx[1]-53.6935750346352*rdx2SqVol[0]*phiC[1]+72.0*rdx2SqVol[0]*bcVals[1]+(3.0*phiLx[0]-39.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]+8.660254037844386*rdx2SqVol[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1]-18.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiLx[1]-131.0*rdx2SqVol[0]*phiC[1]+96.99484522385713*rdx2SqVol[0]*bcVals[1]+((-19.05255888325765*phiLx[0])-29.44486372867091*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+1.732050807568877*rdx2SqVol[0]*phiLx[3]-53.6935750346352*rdx2SqVol[0]*phiC[3]-7.0*rdx2SqVol[1]*phiUy[2]-7.0*rdx2SqVol[1]*phiLy[2]+3.0*rdx2SqVol[0]*phiLx[2]+((-46.0*rdx2SqVol[1])-39.0*rdx2SqVol[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdx2SqVol[1]*phiUy[3]-7.0*rdx2SqVol[1]*phiLy[3]-19.0*rdx2SqVol[0]*phiLx[3]+((-46.0*rdx2SqVol[1])-131.0*rdx2SqVol[0])*phiC[3]-19.05255888325765*rdx2SqVol[0]*phiLx[2]-29.44486372867091*rdx2SqVol[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdx2SqVol[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+8.660254037844386*rdx2SqVol[0]*phiC[1]+16.0*rdx2SqVol[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac-77.94228634059945*rdx2SqVol[1]*phiUy[3]+77.94228634059945*rdx2SqVol[1]*phiLy[3]+(81.0*phiUy[1]+81.0*phiLy[1]-162.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiLx[1]-423.0*rdx2SqVol[0]*phiC[1]+193.9896904477143*rdx2SqVol[0]*bcVals[1]+(98.726896031426*phiC[0]-98.726896031426*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+8.660254037844386*rdx2SqVol[0]*phiLx[3]+8.660254037844386*rdx2SqVol[0]*phiC[3]-7.0*rdx2SqVol[1]*phiUy[2]-7.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-46.0*rdx2SqVol[1])-9.0*rdx2SqVol[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-21.0*rdx2SqVol[1]*phiUy[3]-21.0*rdx2SqVol[1]*phiLy[3]-29.0*rdx2SqVol[0]*phiLx[3]+((-138.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]-32.90896534380867*rdx2SqVol[0]*phiLx[2]+32.90896534380867*rdx2SqVol[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[2]+53.6935750346352*rdx2SqVol[1]*phiC[2]+72.0*rdx2SqVol[1]*bcVals[2]+(3.0*phiUy[0]-39.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[3]+53.6935750346352*rdx2SqVol[1]*phiC[3]+(3.0*phiUy[1]-39.0*phiC[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]+8.660254037844386*rdx2SqVol[0]*phiLx[3]-19.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-131.0*rdx2SqVol[1])-18.0*rdx2SqVol[0])*phiC[2]-96.99484522385713*rdx2SqVol[1]*bcVals[2]+(19.05255888325765*phiUy[0]+29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-19.0*rdx2SqVol[1]*phiUy[3]-7.0*rdx2SqVol[0]*phiUx[3]-7.0*rdx2SqVol[0]*phiLx[3]+((-131.0*rdx2SqVol[1])-46.0*rdx2SqVol[0])*phiC[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+(19.05255888325765*phiUy[1]+29.44486372867091*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]-8.660254037844386*rdx2SqVol[1]*phiC[2]-16.0*rdx2SqVol[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]-8.660254037844386*rdx2SqVol[1]*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac-77.94228634059945*rdx2SqVol[0]*phiUx[3]+77.94228634059945*rdx2SqVol[0]*phiLx[3]-87.0*rdx2SqVol[1]*phiUy[2]+81.0*rdx2SqVol[0]*phiUx[2]+81.0*rdx2SqVol[0]*phiLx[2]+((-423.0*rdx2SqVol[1])-162.0*rdx2SqVol[0])*phiC[2]+193.9896904477143*rdx2SqVol[1]*bcVals[2]+(98.726896031426*phiUy[0]-98.726896031426*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiUy[3]-21.0*rdx2SqVol[0]*phiUx[3]-21.0*rdx2SqVol[0]*phiLx[3]+((-141.0*rdx2SqVol[1])-138.0*rdx2SqVol[0])*phiC[3]+25.98076211353316*rdx2SqVol[0]*phiUx[2]-25.98076211353316*rdx2SqVol[0]*phiLx[2]+(32.90896534380867*phiUy[1]-32.90896534380867*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+72.0*rdx2SqVol[1]*bcVals[3]+1.732050807568877*rdx2SqVol[1]*phiLy[2]-53.6935750346352*rdx2SqVol[1]*phiC[2]+(3.0*phiLy[0]-39.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+1.732050807568877*rdx2SqVol[1]*phiLy[3]-53.6935750346352*rdx2SqVol[1]*phiC[3]+(3.0*phiLy[1]-39.0*phiC[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]+8.660254037844386*rdx2SqVol[0]*phiLx[3]+96.99484522385713*rdx2SqVol[1]*bcVals[3]+9.0*rdx2SqVol[0]*phiUx[2]-19.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-131.0*rdx2SqVol[1])-18.0*rdx2SqVol[0])*phiC[2]+((-19.05255888325765*phiLy[0])-29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdx2SqVol[0]*phiUx[3]-19.0*rdx2SqVol[1]*phiLy[3]-7.0*rdx2SqVol[0]*phiLx[3]+((-131.0*rdx2SqVol[1])-46.0*rdx2SqVol[0])*phiC[3]+8.660254037844386*rdx2SqVol[0]*phiUx[2]-8.660254037844386*rdx2SqVol[0]*phiLx[2]+((-19.05255888325765*phiLy[1])-29.44486372867091*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdx2SqVol[1]*bcVals[3]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+8.660254037844386*rdx2SqVol[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+8.660254037844386*rdx2SqVol[1]*phiLy[3]+8.660254037844386*rdx2SqVol[1]*phiC[3]+(9.0*phiLy[1]-9.0*phiC[1])*rdx2SqVol[1]-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac-77.94228634059945*rdx2SqVol[0]*phiUx[3]+77.94228634059945*rdx2SqVol[0]*phiLx[3]+193.9896904477143*rdx2SqVol[1]*bcVals[3]+81.0*rdx2SqVol[0]*phiUx[2]-87.0*rdx2SqVol[1]*phiLy[2]+81.0*rdx2SqVol[0]*phiLx[2]+((-423.0*rdx2SqVol[1])-162.0*rdx2SqVol[0])*phiC[2]+(98.726896031426*phiC[0]-98.726896031426*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-21.0*rdx2SqVol[0]*phiUx[3]-29.0*rdx2SqVol[1]*phiLy[3]-21.0*rdx2SqVol[0]*phiLx[3]+((-141.0*rdx2SqVol[1])-138.0*rdx2SqVol[0])*phiC[3]+25.98076211353316*rdx2SqVol[0]*phiUx[2]-25.98076211353316*rdx2SqVol[0]*phiLx[2]+(32.90896534380867*phiC[1]-32.90896534380867*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxDirichletLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[2]+53.6935750346352*rdx2SqVol[1]*phiC[2]+72.0*rdx2SqVol[1]*bcVals[2]+(3.0*phiUy[0]-39.0*phiC[0])*rdx2SqVol[1]-1.732050807568877*rdx2SqVol[0]*phiUx[1]+53.6935750346352*rdx2SqVol[0]*phiC[1]+(3.0*phiUx[0]-39.0*phiC[0]+72.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[3]+53.6935750346352*rdx2SqVol[1]*phiC[3]+(3.0*phiUy[1]-39.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiUx[1]-131.0*rdx2SqVol[0]*phiC[1]+(19.05255888325765*phiUx[0]+29.44486372867091*phiC[0]-96.99484522385713*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-1.732050807568877*rdx2SqVol[0]*phiUx[3]+53.6935750346352*rdx2SqVol[0]*phiC[3]-19.0*rdx2SqVol[1]*phiUy[2]+3.0*rdx2SqVol[0]*phiUx[2]+((-131.0*rdx2SqVol[1])-39.0*rdx2SqVol[0])*phiC[2]-96.99484522385713*rdx2SqVol[1]*bcVals[2]+(19.05255888325765*phiUy[0]+29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-19.0*rdx2SqVol[1]*phiUy[3]-19.0*rdx2SqVol[0]*phiUx[3]+((-131.0*rdx2SqVol[1])-131.0*rdx2SqVol[0])*phiC[3]+19.05255888325765*rdx2SqVol[0]*phiUx[2]+29.44486372867091*rdx2SqVol[0]*phiC[2]+(19.05255888325765*phiUy[1]+29.44486372867091*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxDirichletLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]-8.660254037844386*rdx2SqVol[1]*phiC[2]-16.0*rdx2SqVol[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdx2SqVol[1]-1.732050807568877*rdx2SqVol[0]*phiUx[1]+53.6935750346352*rdx2SqVol[0]*phiC[1]+(3.0*phiUx[0]-39.0*phiC[0]+72.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]-8.660254037844386*rdx2SqVol[1]*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiUx[1]-131.0*rdx2SqVol[0]*phiC[1]+(19.05255888325765*phiUx[0]+29.44486372867091*phiC[0]-96.99484522385713*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac-15.58845726811989*rdx2SqVol[0]*phiUx[3]+483.2421753117166*rdx2SqVol[0]*phiC[3]-87.0*rdx2SqVol[1]*phiUy[2]+27.0*rdx2SqVol[0]*phiUx[2]+((-423.0*rdx2SqVol[1])-351.0*rdx2SqVol[0])*phiC[2]+193.9896904477143*rdx2SqVol[1]*bcVals[2]+(98.726896031426*phiUy[0]-98.726896031426*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiUy[3]-57.0*rdx2SqVol[0]*phiUx[3]+((-141.0*rdx2SqVol[1])-393.0*rdx2SqVol[0])*phiC[3]+57.15767664977295*rdx2SqVol[0]*phiUx[2]+88.33459118601273*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiUy[1]-32.90896534380867*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxNeumannLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[2]+53.6935750346352*rdx2SqVol[1]*phiC[2]+72.0*rdx2SqVol[1]*bcVals[2]+(3.0*phiUy[0]-39.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]-8.660254037844386*rdx2SqVol[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac-15.58845726811989*rdx2SqVol[1]*phiUy[3]+483.2421753117166*rdx2SqVol[1]*phiC[3]+(27.0*phiUy[1]-351.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiUx[1]-423.0*rdx2SqVol[0]*phiC[1]+(98.726896031426*phiUx[0]-98.726896031426*phiC[0]+193.9896904477143*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]-8.660254037844386*rdx2SqVol[0]*phiC[3]-19.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiUx[2]+((-131.0*rdx2SqVol[1])-9.0*rdx2SqVol[0])*phiC[2]-96.99484522385713*rdx2SqVol[1]*bcVals[2]+(19.05255888325765*phiUy[0]+29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-57.0*rdx2SqVol[1]*phiUy[3]-29.0*rdx2SqVol[0]*phiUx[3]+((-393.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]+32.90896534380867*rdx2SqVol[0]*phiUx[2]-32.90896534380867*rdx2SqVol[0]*phiC[2]+(57.15767664977295*phiUy[1]+88.33459118601273*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxNeumannLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]-8.660254037844386*rdx2SqVol[1]*phiC[2]-16.0*rdx2SqVol[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]-8.660254037844386*rdx2SqVol[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac-77.94228634059945*rdx2SqVol[1]*phiUy[3]-77.94228634059945*rdx2SqVol[1]*phiC[3]+(81.0*phiUy[1]-81.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiUx[1]-423.0*rdx2SqVol[0]*phiC[1]+(98.726896031426*phiUx[0]-98.726896031426*phiC[0]+193.9896904477143*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac-77.94228634059945*rdx2SqVol[0]*phiUx[3]-77.94228634059945*rdx2SqVol[0]*phiC[3]-87.0*rdx2SqVol[1]*phiUy[2]+81.0*rdx2SqVol[0]*phiUx[2]+((-423.0*rdx2SqVol[1])-81.0*rdx2SqVol[0])*phiC[2]+193.9896904477143*rdx2SqVol[1]*bcVals[2]+(98.726896031426*phiUy[0]-98.726896031426*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiUy[3]-29.0*rdx2SqVol[0]*phiUx[3]+((-141.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]+32.90896534380867*rdx2SqVol[0]*phiUx[2]-32.90896534380867*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiUy[1]-32.90896534380867*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxDirichletUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+72.0*rdx2SqVol[1]*bcVals[3]+1.732050807568877*rdx2SqVol[1]*phiLy[2]-53.6935750346352*rdx2SqVol[1]*phiC[2]+(3.0*phiLy[0]-39.0*phiC[0])*rdx2SqVol[1]-1.732050807568877*rdx2SqVol[0]*phiUx[1]+53.6935750346352*rdx2SqVol[0]*phiC[1]+(3.0*phiUx[0]-39.0*phiC[0]+72.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+1.732050807568877*rdx2SqVol[1]*phiLy[3]-53.6935750346352*rdx2SqVol[1]*phiC[3]+(3.0*phiLy[1]-39.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiUx[1]-131.0*rdx2SqVol[0]*phiC[1]+(19.05255888325765*phiUx[0]+29.44486372867091*phiC[0]-96.99484522385713*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-1.732050807568877*rdx2SqVol[0]*phiUx[3]+53.6935750346352*rdx2SqVol[0]*phiC[3]+96.99484522385713*rdx2SqVol[1]*bcVals[3]+3.0*rdx2SqVol[0]*phiUx[2]-19.0*rdx2SqVol[1]*phiLy[2]+((-131.0*rdx2SqVol[1])-39.0*rdx2SqVol[0])*phiC[2]+((-19.05255888325765*phiLy[0])-29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-19.0*rdx2SqVol[0]*phiUx[3]-19.0*rdx2SqVol[1]*phiLy[3]+((-131.0*rdx2SqVol[1])-131.0*rdx2SqVol[0])*phiC[3]+19.05255888325765*rdx2SqVol[0]*phiUx[2]+29.44486372867091*rdx2SqVol[0]*phiC[2]+((-19.05255888325765*phiLy[1])-29.44486372867091*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxDirichletUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdx2SqVol[1]*bcVals[3]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+8.660254037844386*rdx2SqVol[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdx2SqVol[1]-1.732050807568877*rdx2SqVol[0]*phiUx[1]+53.6935750346352*rdx2SqVol[0]*phiC[1]+(3.0*phiUx[0]-39.0*phiC[0]+72.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+8.660254037844386*rdx2SqVol[1]*phiLy[3]+8.660254037844386*rdx2SqVol[1]*phiC[3]+(9.0*phiLy[1]-9.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiUx[1]-131.0*rdx2SqVol[0]*phiC[1]+(19.05255888325765*phiUx[0]+29.44486372867091*phiC[0]-96.99484522385713*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac-15.58845726811989*rdx2SqVol[0]*phiUx[3]+483.2421753117166*rdx2SqVol[0]*phiC[3]+193.9896904477143*rdx2SqVol[1]*bcVals[3]+27.0*rdx2SqVol[0]*phiUx[2]-87.0*rdx2SqVol[1]*phiLy[2]+((-423.0*rdx2SqVol[1])-351.0*rdx2SqVol[0])*phiC[2]+(98.726896031426*phiC[0]-98.726896031426*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-57.0*rdx2SqVol[0]*phiUx[3]-29.0*rdx2SqVol[1]*phiLy[3]+((-141.0*rdx2SqVol[1])-393.0*rdx2SqVol[0])*phiC[3]+57.15767664977295*rdx2SqVol[0]*phiUx[2]+88.33459118601273*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiC[1]-32.90896534380867*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxNeumannUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+72.0*rdx2SqVol[1]*bcVals[3]+1.732050807568877*rdx2SqVol[1]*phiLy[2]-53.6935750346352*rdx2SqVol[1]*phiC[2]+(3.0*phiLy[0]-39.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]-8.660254037844386*rdx2SqVol[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac+15.58845726811989*rdx2SqVol[1]*phiLy[3]-483.2421753117166*rdx2SqVol[1]*phiC[3]+(27.0*phiLy[1]-351.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiUx[1]-423.0*rdx2SqVol[0]*phiC[1]+(98.726896031426*phiUx[0]-98.726896031426*phiC[0]+193.9896904477143*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[3]-8.660254037844386*rdx2SqVol[0]*phiC[3]+96.99484522385713*rdx2SqVol[1]*bcVals[3]+9.0*rdx2SqVol[0]*phiUx[2]-19.0*rdx2SqVol[1]*phiLy[2]+((-131.0*rdx2SqVol[1])-9.0*rdx2SqVol[0])*phiC[2]+((-19.05255888325765*phiLy[0])-29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[0]*phiUx[3]-57.0*rdx2SqVol[1]*phiLy[3]+((-393.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]+32.90896534380867*rdx2SqVol[0]*phiUx[2]-32.90896534380867*rdx2SqVol[0]*phiC[2]+((-57.15767664977295*phiLy[1])-88.33459118601273*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_LxNeumannUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdx2SqVol[1]*bcVals[3]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+8.660254037844386*rdx2SqVol[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdx2SqVol[1]-8.660254037844386*rdx2SqVol[0]*phiUx[1]-8.660254037844386*rdx2SqVol[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac+77.94228634059945*rdx2SqVol[1]*phiLy[3]+77.94228634059945*rdx2SqVol[1]*phiC[3]+(81.0*phiLy[1]-81.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiUx[1]-423.0*rdx2SqVol[0]*phiC[1]+(98.726896031426*phiUx[0]-98.726896031426*phiC[0]+193.9896904477143*bcVals[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac-77.94228634059945*rdx2SqVol[0]*phiUx[3]-77.94228634059945*rdx2SqVol[0]*phiC[3]+193.9896904477143*rdx2SqVol[1]*bcVals[3]+81.0*rdx2SqVol[0]*phiUx[2]-87.0*rdx2SqVol[1]*phiLy[2]+((-423.0*rdx2SqVol[1])-81.0*rdx2SqVol[0])*phiC[2]+(98.726896031426*phiC[0]-98.726896031426*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[0]*phiUx[3]-29.0*rdx2SqVol[1]*phiLy[3]+((-141.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]+32.90896534380867*rdx2SqVol[0]*phiUx[2]-32.90896534380867*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiC[1]-32.90896534380867*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxDirichletLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[2]+53.6935750346352*rdx2SqVol[1]*phiC[2]+72.0*rdx2SqVol[1]*bcVals[2]+(3.0*phiUy[0]-39.0*phiC[0])*rdx2SqVol[1]+1.732050807568877*rdx2SqVol[0]*phiLx[1]-53.6935750346352*rdx2SqVol[0]*phiC[1]+72.0*rdx2SqVol[0]*bcVals[1]+(3.0*phiLx[0]-39.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[3]+53.6935750346352*rdx2SqVol[1]*phiC[3]+(3.0*phiUy[1]-39.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiLx[1]-131.0*rdx2SqVol[0]*phiC[1]+96.99484522385713*rdx2SqVol[0]*bcVals[1]+((-19.05255888325765*phiLx[0])-29.44486372867091*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+1.732050807568877*rdx2SqVol[0]*phiLx[3]-53.6935750346352*rdx2SqVol[0]*phiC[3]-19.0*rdx2SqVol[1]*phiUy[2]+3.0*rdx2SqVol[0]*phiLx[2]+((-131.0*rdx2SqVol[1])-39.0*rdx2SqVol[0])*phiC[2]-96.99484522385713*rdx2SqVol[1]*bcVals[2]+(19.05255888325765*phiUy[0]+29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-19.0*rdx2SqVol[1]*phiUy[3]-19.0*rdx2SqVol[0]*phiLx[3]+((-131.0*rdx2SqVol[1])-131.0*rdx2SqVol[0])*phiC[3]-19.05255888325765*rdx2SqVol[0]*phiLx[2]-29.44486372867091*rdx2SqVol[0]*phiC[2]+(19.05255888325765*phiUy[1]+29.44486372867091*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxDirichletLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]-8.660254037844386*rdx2SqVol[1]*phiC[2]-16.0*rdx2SqVol[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdx2SqVol[1]+1.732050807568877*rdx2SqVol[0]*phiLx[1]-53.6935750346352*rdx2SqVol[0]*phiC[1]+72.0*rdx2SqVol[0]*bcVals[1]+(3.0*phiLx[0]-39.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[3]-8.660254037844386*rdx2SqVol[1]*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiLx[1]-131.0*rdx2SqVol[0]*phiC[1]+96.99484522385713*rdx2SqVol[0]*bcVals[1]+((-19.05255888325765*phiLx[0])-29.44486372867091*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac+15.58845726811989*rdx2SqVol[0]*phiLx[3]-483.2421753117166*rdx2SqVol[0]*phiC[3]-87.0*rdx2SqVol[1]*phiUy[2]+27.0*rdx2SqVol[0]*phiLx[2]+((-423.0*rdx2SqVol[1])-351.0*rdx2SqVol[0])*phiC[2]+193.9896904477143*rdx2SqVol[1]*bcVals[2]+(98.726896031426*phiUy[0]-98.726896031426*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiUy[3]-57.0*rdx2SqVol[0]*phiLx[3]+((-141.0*rdx2SqVol[1])-393.0*rdx2SqVol[0])*phiC[3]-57.15767664977295*rdx2SqVol[0]*phiLx[2]-88.33459118601273*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiUy[1]-32.90896534380867*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxNeumannLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-1.732050807568877*rdx2SqVol[1]*phiUy[2]+53.6935750346352*rdx2SqVol[1]*phiC[2]+72.0*rdx2SqVol[1]*bcVals[2]+(3.0*phiUy[0]-39.0*phiC[0])*rdx2SqVol[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+8.660254037844386*rdx2SqVol[0]*phiC[1]+16.0*rdx2SqVol[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac-15.58845726811989*rdx2SqVol[1]*phiUy[3]+483.2421753117166*rdx2SqVol[1]*phiC[3]+(27.0*phiUy[1]-351.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiLx[1]-423.0*rdx2SqVol[0]*phiC[1]+193.9896904477143*rdx2SqVol[0]*bcVals[1]+(98.726896031426*phiC[0]-98.726896031426*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+8.660254037844386*rdx2SqVol[0]*phiLx[3]+8.660254037844386*rdx2SqVol[0]*phiC[3]-19.0*rdx2SqVol[1]*phiUy[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-131.0*rdx2SqVol[1])-9.0*rdx2SqVol[0])*phiC[2]-96.99484522385713*rdx2SqVol[1]*bcVals[2]+(19.05255888325765*phiUy[0]+29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-57.0*rdx2SqVol[1]*phiUy[3]-29.0*rdx2SqVol[0]*phiLx[3]+((-393.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]-32.90896534380867*rdx2SqVol[0]*phiLx[2]+32.90896534380867*rdx2SqVol[0]*phiC[2]+(57.15767664977295*phiUy[1]+88.33459118601273*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxNeumannLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[1]*phiUy[2]-8.660254037844386*rdx2SqVol[1]*phiC[2]-16.0*rdx2SqVol[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdx2SqVol[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+8.660254037844386*rdx2SqVol[0]*phiC[1]+16.0*rdx2SqVol[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac-77.94228634059945*rdx2SqVol[1]*phiUy[3]-77.94228634059945*rdx2SqVol[1]*phiC[3]+(81.0*phiUy[1]-81.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiLx[1]-423.0*rdx2SqVol[0]*phiC[1]+193.9896904477143*rdx2SqVol[0]*bcVals[1]+(98.726896031426*phiC[0]-98.726896031426*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac+77.94228634059945*rdx2SqVol[0]*phiLx[3]+77.94228634059945*rdx2SqVol[0]*phiC[3]-87.0*rdx2SqVol[1]*phiUy[2]+81.0*rdx2SqVol[0]*phiLx[2]+((-423.0*rdx2SqVol[1])-81.0*rdx2SqVol[0])*phiC[2]+193.9896904477143*rdx2SqVol[1]*bcVals[2]+(98.726896031426*phiUy[0]-98.726896031426*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiUy[3]-29.0*rdx2SqVol[0]*phiLx[3]+((-141.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]-32.90896534380867*rdx2SqVol[0]*phiLx[2]+32.90896534380867*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiUy[1]-32.90896534380867*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxDirichletUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+72.0*rdx2SqVol[1]*bcVals[3]+1.732050807568877*rdx2SqVol[1]*phiLy[2]-53.6935750346352*rdx2SqVol[1]*phiC[2]+(3.0*phiLy[0]-39.0*phiC[0])*rdx2SqVol[1]+1.732050807568877*rdx2SqVol[0]*phiLx[1]-53.6935750346352*rdx2SqVol[0]*phiC[1]+72.0*rdx2SqVol[0]*bcVals[1]+(3.0*phiLx[0]-39.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+1.732050807568877*rdx2SqVol[1]*phiLy[3]-53.6935750346352*rdx2SqVol[1]*phiC[3]+(3.0*phiLy[1]-39.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiLx[1]-131.0*rdx2SqVol[0]*phiC[1]+96.99484522385713*rdx2SqVol[0]*bcVals[1]+((-19.05255888325765*phiLx[0])-29.44486372867091*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+1.732050807568877*rdx2SqVol[0]*phiLx[3]-53.6935750346352*rdx2SqVol[0]*phiC[3]+96.99484522385713*rdx2SqVol[1]*bcVals[3]-19.0*rdx2SqVol[1]*phiLy[2]+3.0*rdx2SqVol[0]*phiLx[2]+((-131.0*rdx2SqVol[1])-39.0*rdx2SqVol[0])*phiC[2]+((-19.05255888325765*phiLy[0])-29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-19.0*rdx2SqVol[1]*phiLy[3]-19.0*rdx2SqVol[0]*phiLx[3]+((-131.0*rdx2SqVol[1])-131.0*rdx2SqVol[0])*phiC[3]-19.05255888325765*rdx2SqVol[0]*phiLx[2]-29.44486372867091*rdx2SqVol[0]*phiC[2]+((-19.05255888325765*phiLy[1])-29.44486372867091*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxDirichletUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdx2SqVol[1]*bcVals[3]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+8.660254037844386*rdx2SqVol[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdx2SqVol[1]+1.732050807568877*rdx2SqVol[0]*phiLx[1]-53.6935750346352*rdx2SqVol[0]*phiC[1]+72.0*rdx2SqVol[0]*bcVals[1]+(3.0*phiLx[0]-39.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+8.660254037844386*rdx2SqVol[1]*phiLy[3]+8.660254037844386*rdx2SqVol[1]*phiC[3]+(9.0*phiLy[1]-9.0*phiC[1])*rdx2SqVol[1]-19.0*rdx2SqVol[0]*phiLx[1]-131.0*rdx2SqVol[0]*phiC[1]+96.99484522385713*rdx2SqVol[0]*bcVals[1]+((-19.05255888325765*phiLx[0])-29.44486372867091*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac+15.58845726811989*rdx2SqVol[0]*phiLx[3]-483.2421753117166*rdx2SqVol[0]*phiC[3]+193.9896904477143*rdx2SqVol[1]*bcVals[3]-87.0*rdx2SqVol[1]*phiLy[2]+27.0*rdx2SqVol[0]*phiLx[2]+((-423.0*rdx2SqVol[1])-351.0*rdx2SqVol[0])*phiC[2]+(98.726896031426*phiC[0]-98.726896031426*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiLy[3]-57.0*rdx2SqVol[0]*phiLx[3]+((-141.0*rdx2SqVol[1])-393.0*rdx2SqVol[0])*phiC[3]-57.15767664977295*rdx2SqVol[0]*phiLx[2]-88.33459118601273*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiC[1]-32.90896534380867*phiLy[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxNeumannUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+72.0*rdx2SqVol[1]*bcVals[3]+1.732050807568877*rdx2SqVol[1]*phiLy[2]-53.6935750346352*rdx2SqVol[1]*phiC[2]+(3.0*phiLy[0]-39.0*phiC[0])*rdx2SqVol[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+8.660254037844386*rdx2SqVol[0]*phiC[1]+16.0*rdx2SqVol[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac+15.58845726811989*rdx2SqVol[1]*phiLy[3]-483.2421753117166*rdx2SqVol[1]*phiC[3]+(27.0*phiLy[1]-351.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiLx[1]-423.0*rdx2SqVol[0]*phiC[1]+193.9896904477143*rdx2SqVol[0]*bcVals[1]+(98.726896031426*phiC[0]-98.726896031426*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+8.660254037844386*rdx2SqVol[0]*phiLx[3]+8.660254037844386*rdx2SqVol[0]*phiC[3]+96.99484522385713*rdx2SqVol[1]*bcVals[3]-19.0*rdx2SqVol[1]*phiLy[2]+9.0*rdx2SqVol[0]*phiLx[2]+((-131.0*rdx2SqVol[1])-9.0*rdx2SqVol[0])*phiC[2]+((-19.05255888325765*phiLy[0])-29.44486372867091*phiC[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-57.0*rdx2SqVol[1]*phiLy[3]-29.0*rdx2SqVol[0]*phiLx[3]+((-393.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]-32.90896534380867*rdx2SqVol[0]*phiLx[2]+32.90896534380867*rdx2SqVol[0]*phiC[2]+((-57.15767664977295*phiLy[1])-88.33459118601273*phiC[1])*rdx2SqVol[1]); 

}

void MGpoissonResidue2xSer_UxNeumannUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  double rdx2SqVolSq[2]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVolSq[0] = rdx2SqVol[0]*rdx2SqVol[0]; 
  rdx2SqVol[1]   = volFac*4.0/(dxC[1]*dxC[1]); 
  rdx2SqVolSq[1] = rdx2SqVol[1]*rdx2SqVol[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdx2SqVol[1]*bcVals[3]+8.660254037844386*rdx2SqVol[1]*phiLy[2]+8.660254037844386*rdx2SqVol[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdx2SqVol[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+8.660254037844386*rdx2SqVol[0]*phiC[1]+16.0*rdx2SqVol[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.006944444444444444*(144.0*rho[1]*volFac+77.94228634059945*rdx2SqVol[1]*phiLy[3]+77.94228634059945*rdx2SqVol[1]*phiC[3]+(81.0*phiLy[1]-81.0*phiC[1])*rdx2SqVol[1]-87.0*rdx2SqVol[0]*phiLx[1]-423.0*rdx2SqVol[0]*phiC[1]+193.9896904477143*rdx2SqVol[0]*bcVals[1]+(98.726896031426*phiC[0]-98.726896031426*phiLx[0])*rdx2SqVol[0]); 
  resOut[2] = 0.006944444444444444*(144.0*rho[2]*volFac+77.94228634059945*rdx2SqVol[0]*phiLx[3]+77.94228634059945*rdx2SqVol[0]*phiC[3]+193.9896904477143*rdx2SqVol[1]*bcVals[3]-87.0*rdx2SqVol[1]*phiLy[2]+81.0*rdx2SqVol[0]*phiLx[2]+((-423.0*rdx2SqVol[1])-81.0*rdx2SqVol[0])*phiC[2]+(98.726896031426*phiC[0]-98.726896031426*phiLy[0])*rdx2SqVol[1]); 
  resOut[3] = 0.02083333333333333*(48.0*rho[3]*volFac-29.0*rdx2SqVol[1]*phiLy[3]-29.0*rdx2SqVol[0]*phiLx[3]+((-141.0*rdx2SqVol[1])-141.0*rdx2SqVol[0])*phiC[3]-32.90896534380867*rdx2SqVol[0]*phiLx[2]+32.90896534380867*rdx2SqVol[0]*phiC[2]+(32.90896534380867*phiC[1]-32.90896534380867*phiLy[1])*rdx2SqVol[1]); 

}

