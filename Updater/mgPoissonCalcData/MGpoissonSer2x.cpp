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
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 
  double *dxUy = dx[3]; 
  double *dxLy = dx[4]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdxCp2[2]; 
  double rdxLx[2]; 
  double rdxUx[2]; 
  double rdxLxSq[2]; 
  double rdxUxSq[2]; 
  double rdxLxCu[2]; 
  double rdxUxCu[2]; 
  double rdxLxR4[2]; 
  double rdxUxR4[2]; 
  double rdxLy[2]; 
  double rdxUy[2]; 
  double rdxLySq[2]; 
  double rdxUySq[2]; 
  double rdxLyCu[2]; 
  double rdxUyCu[2]; 
  double rdxLyR4[2]; 
  double rdxUyR4[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = volFac*4.0/(dxLx[1]*dxLx[1]); 
  rdxUx[1]   = volFac*4.0/(dxUx[1]*dxUx[1]); 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = volFac*4.0/(dxLy[0]*dxLy[0]); 
  rdxUy[0]   = volFac*4.0/(dxUy[0]*dxUy[0]); 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = volFac*4.0/(dxLy[1]*dxLy[1]); 
  rdxUy[1]   = volFac*4.0/(dxUy[1]*dxUy[1]); 
  rdxLySq[1] = rdxLy[1]*rdxLy[1]; 
  rdxUySq[1] = rdxUy[1]*rdxUy[1]; 
  rdxLyCu[1] = rdxLySq[1]*rdxLy[1]; 
  rdxUyCu[1] = rdxUySq[1]*rdxUy[1]; 
  rdxLyR4[1] = rdxLyCu[1]*rdxLy[1]; 
  rdxUyR4[1] = rdxUyCu[1]*rdxUy[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((((9600.0*rdxUx[0]-9600.0*rdxLx[0])*rdxUySq[1]+(9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUy[1]+(9600.0*rdxLx[0]-9600.0*rdxUx[0])*rdxLySq[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLy[1])*rho[3]+((-415.6921938165305*rdxUyCu[1])+((-27435.68479189101*rdxLy[1])-25495.78788741387*rdxUx[0]-25495.78788741387*rdxLx[0])*rdxUySq[1]+(27435.68479189101*rdxLySq[1]-25080.09569359734*rdxUxSq[0]-23140.1987891202*rdxLx[0]*rdxUx[0]-25080.09569359734*rdxLxSq[0])*rdxUy[1]+415.6921938165305*rdxLyCu[1]+(25495.78788741387*rdxUx[0]+25495.78788741387*rdxLx[0])*rdxLySq[1]+(25080.09569359734*rdxUxSq[0]+23140.1987891202*rdxLx[0]*rdxUx[0]+25080.09569359734*rdxLxSq[0])*rdxLy[1])*rho[2]+((25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxUySq[1]+((23140.1987891202*rdxLx[0]-23140.1987891202*rdxUx[0])*rdxLy[1]-25495.78788741387*rdxUxSq[0]+25495.78788741387*rdxLxSq[0])*rdxUy[1]+(25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxLySq[1]+(25495.78788741387*rdxLxSq[0]-25495.78788741387*rdxUxSq[0])*rdxLy[1]-415.6921938165305*rdxUxCu[0]-27435.68479189101*rdxLx[0]*rdxUxSq[0]+27435.68479189101*rdxLxSq[0]*rdxUx[0]+415.6921938165305*rdxLxCu[0])*rho[1]+1104.0*rho[0]*rdxUyCu[1]+(75072.0*rho[0]*rdxLy[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0])*rdxUySq[1]+(75072.0*rho[0]*rdxLySq[1]+(164368.0*rdxUx[0]+164368.0*rdxLx[0])*rho[0]*rdxLy[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0])*rdxUy[1]+1104.0*rho[0]*rdxLyCu[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0]*rdxLySq[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0]*rdxLy[1]+(1104.0*rdxUxCu[0]+75072.0*rdxLx[0]*rdxUxSq[0]+75072.0*rdxLxSq[0]*rdxUx[0]+1104.0*rdxLxCu[0])*rho[0])*volFac+((9375.0*rdxUx[0]-9375.0*rdxLx[0])*rdxUyCu[1]+((12525.0*rdxUx[0]-12525.0*rdxLx[0])*rdxLy[1]+9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+((17775.0*rdxUx[0]-17775.0*rdxLx[0])*rdxLySq[1]+(18000.0*rdxUxSq[0]-18000.0*rdxLxSq[0])*rdxLy[1]+225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]+9600.0*rdxUxSq[0]+18000.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])+9375.0*rdxUxCu[0]+12525.0*rdxLx[0]*rdxUxSq[0]+17775.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+((-9600.0*rdxUxSq[0])-18000.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-9375.0*rdxUxCu[0])-12525.0*rdxLx[0]*rdxUxSq[0]-17775.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((17775.0*rdxLx[0]-17775.0*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((12525.0*rdxLx[0]-12525.0*rdxUx[0])*rdxLySq[1]+(18000.0*rdxLxSq[0]-18000.0*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9375.0*rdxLx[0]-9375.0*rdxUx[0])*rdxLyCu[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLySq[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-18000.0*rdxLx[0]*rdxUx[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-17775.0*rdxLx[0]*rdxUxSq[0]-12525.0*rdxLxSq[0]*rdxUx[0]-9375.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(18000.0*rdxLx[0]*rdxUx[0]+9600.0*rdxLxSq[0])*rdxLySq[1]+(17775.0*rdxLx[0]*rdxUxSq[0]+12525.0*rdxLxSq[0]*rdxUx[0]+9375.0*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((-415.6921938165305*rdxUyR4[1])+((-28630.79984911354*rdxLy[1])-25729.61474643567*rdxUx[0]-25729.61474643567*rdxLx[0])*rdxUyCu[1]+((-52637.02404201817*rdxLySq[1])+((-88966.78973077537*rdxUx[0])-88966.78973077537*rdxLx[0])*rdxLy[1]-25911.4800812304*rdxUxSq[0]-78842.95276053529*rdxLx[0]*rdxUx[0]-25911.4800812304*rdxLxSq[0])*rdxUySq[1]+((-779.4228634059946*rdxLyCu[1])+((-48038.4291479228*rdxUx[0])-48038.4291479228*rdxLx[0])*rdxLySq[1]+((-47856.56381312806*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0]-47856.56381312806*rdxLxSq[0])*rdxLy[1]-597.5575286112626*rdxUxCu[0]-40633.91194556586*rdxLx[0]*rdxUxSq[0]-40633.91194556586*rdxLxSq[0]*rdxUx[0]-597.5575286112626*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-233.8268590217983*rdxUx[0]*rdxUyCu[1])+((-15432.57269543869*rdxUx[0]*rdxLy[1])-9145.22826396367*rdxUxSq[0]-19537.53310937693*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15432.57269543869*rdxUx[0]*rdxLySq[1]-8911.401404941873*rdxUxCu[0]-13016.36181888011*rdxLx[0]*rdxUxSq[0]-19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]+19537.53310937693*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8911.401404941873*rdxUxCu[0]+13016.36181888011*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+(779.4228634059946*rdxLy[1]*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(48038.4291479228*rdxUx[0]+48038.4291479228*rdxLx[0])*rdxLy[1])*rdxUySq[1]+(28630.79984911354*rdxLyCu[1]+(88966.78973077537*rdxUx[0]+88966.78973077537*rdxLx[0])*rdxLySq[1]+(47856.56381312806*rdxUxSq[0]+99090.62670101546*rdxLx[0]*rdxUx[0]+47856.56381312806*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+415.6921938165305*rdxLyR4[1]+(25729.61474643567*rdxUx[0]+25729.61474643567*rdxLx[0])*rdxLyCu[1]+(25911.4800812304*rdxUxSq[0]+78842.95276053529*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*rdxLySq[1]+(597.5575286112626*rdxUxCu[0]+40633.91194556586*rdxLx[0]*rdxUxSq[0]+40633.91194556586*rdxLxSq[0]*rdxUx[0]+597.5575286112626*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15432.57269543869*rdxLx[0]*rdxLy[1])-19537.53310937693*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(15432.57269543869*rdxLx[0]*rdxLySq[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-13016.36181888011*rdxLxSq[0]*rdxUx[0]-8911.401404941873*rdxLxCu[0])*rdxUy[1]+233.8268590217983*rdxLx[0]*rdxLyCu[1]+(19537.53310937693*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(19303.70625035514*rdxLx[0]*rdxUxSq[0]+13016.36181888011*rdxLxSq[0]*rdxUx[0]+8911.401404941873*rdxLxCu[0])*rdxLy[1])*phiLx[2]+396.0*phiUy[0]*rdxUyR4[1]+((27378.0*phiUy[0]+846.0*phiLy[0])*rdxLy[1]+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiUy[1]-597.5575286112626*rdxUx[0]*phiUx[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(24531.0*phiUy[0]+621.0*phiUx[0])*rdxUx[0]+(24531.0*phiUy[0]+621.0*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((57078.0*phiUy[0]+57078.0*phiLy[0])*rdxLySq[1]+((13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(92457.0*phiUy[0]+42228.0*phiUx[0]+52131.0*phiLy[0])*rdxUx[0]+(92457.0*phiUy[0]+52131.0*phiLy[0]+42228.0*phiLx[0])*rdxLx[0])*rdxLy[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiUy[1]+((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUy[0]+24756.0*phiUx[0])*rdxUxSq[0]+(79932.0*phiUy[0]+51906.0*phiUx[0]+51906.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiUy[0]+24756.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((846.0*phiUy[0]+27378.0*phiLy[0])*rdxLyCu[1]+((19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(52131.0*phiUy[0]+42228.0*phiUx[0]+92457.0*phiLy[0])*rdxUx[0]+(52131.0*phiUy[0]+92457.0*phiLy[0]+42228.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiUy[1]+((-78842.95276053529*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0])*phiUx[1]+(19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiLy[1]+(99090.62670101546*rdxLx[0]*rdxUx[0]+78842.95276053529*rdxLxSq[0])*phiLx[1]+(51906.0*phiUy[0]+79932.0*phiUx[0]+51906.0*phiLy[0])*rdxUxSq[0]+(104982.0*phiUy[0]+104982.0*phiUx[0]+104982.0*phiLy[0]+104982.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(51906.0*phiUy[0]+51906.0*phiLy[0]+79932.0*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiUy[1]+((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(621.0*phiUy[0]+24531.0*phiUx[0])*rdxUxCu[0]+(42228.0*phiUy[0]+92457.0*phiUx[0]+52131.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(42228.0*phiUy[0]+52131.0*phiUx[0]+92457.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiUy[0]+24531.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+396.0*phiLy[0]*rdxLyR4[1]+((-597.5575286112626*rdxUx[0]*phiUx[1])+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiLy[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(621.0*phiUx[0]+24531.0*phiLy[0])*rdxUx[0]+(24531.0*phiLy[0]+621.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+(((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiLy[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUx[0]+24756.0*phiLy[0])*rdxUxSq[0]+(51906.0*phiUx[0]+79932.0*phiLy[0]+51906.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiLy[0]+24756.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiLy[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(24531.0*phiUx[0]+621.0*phiLy[0])*rdxUxCu[0]+(92457.0*phiUx[0]+42228.0*phiLy[0]+52131.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(52131.0*phiUx[0]+42228.0*phiLy[0]+92457.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiLy[0]+24531.0*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-415.6921938165305*rdxUxR4[0])-28630.79984911354*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(779.4228634059946*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+28630.79984911354*rdxLxCu[0]*rdxUx[0]+415.6921938165305*rdxLxR4[0])*phiLx[1]+396.0*phiUx[0]*rdxUxR4[0]+(27378.0*phiUx[0]+846.0*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(57078.0*phiUx[0]+57078.0*phiLx[0])*rdxLxSq[0]*rdxUxSq[0]+(846.0*phiUx[0]+27378.0*phiLx[0])*rdxLxCu[0]*rdxUx[0]+396.0*phiLx[0]*rdxLxR4[0])/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxUyCu[1]+(27435.68479189101*rdxLy[1]+9976.61265159673*rdxUx[0]+9976.61265159673*rdxLx[0])*rdxUySq[1]+((-27435.68479189101*rdxLySq[1])+9560.920457780201*rdxUxSq[0]-7898.15168251408*rdxLx[0]*rdxUx[0]+9560.920457780201*rdxLxSq[0])*rdxUy[1]-415.6921938165305*rdxLyCu[1]+((-9976.61265159673*rdxUx[0])-9976.61265159673*rdxLx[0])*rdxLySq[1]+((-9560.920457780201*rdxUxSq[0])+7898.15168251408*rdxLx[0]*rdxUx[0]-9560.920457780201*rdxLxSq[0])*rdxLy[1])*rho[3]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[2]+((-1104.0*rdxUyCu[1])+((-75072.0*rdxLy[1])-27600.0*rdxUx[0]-27600.0*rdxLx[0])*rdxUySq[1]+((-75072.0*rdxLySq[1])+((-126960.0*rdxUx[0])-126960.0*rdxLx[0])*rdxLy[1]-26928.0*rdxUxSq[0]-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxUy[1]-1104.0*rdxLyCu[1]+((-27600.0*rdxUx[0])-27600.0*rdxLx[0])*rdxLySq[1]+((-26928.0*rdxUxSq[0])-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxLy[1]-432.0*rdxUxCu[0]-29376.0*rdxLx[0]*rdxUxSq[0]-29376.0*rdxLxSq[0]*rdxUx[0]-432.0*rdxLxCu[0])*rho[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxUySq[1]+((60164.51685171252*rdxUx[0]-60164.51685171252*rdxLx[0])*rho[0]*rdxLy[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0])*rdxUy[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxLySq[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0]*rdxLy[1]+(1080.799703922979*rdxUxCu[0]+71332.78045891662*rdxLx[0]*rdxUxSq[0]-71332.78045891662*rdxLxSq[0]*rdxUx[0]-1080.799703922979*rdxLxCu[0])*rho[0])*volFac+(415.6921938165305*rdxUyR4[1]+(28630.79984911354*rdxLy[1]+10574.170180208*rdxUx[0]+10574.170180208*rdxLx[0])*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(68719.1157902952*rdxUx[0]+68719.1157902952*rdxLx[0])*rdxLy[1]+10392.30484541326*rdxUxSq[0]+47804.60228890101*rdxLx[0]*rdxUx[0]+10392.30484541326*rdxLxSq[0])*rdxUySq[1]+(779.4228634059946*rdxLyCu[1]+(19303.70625035514*rdxUx[0]+19303.70625035514*rdxLx[0])*rdxLySq[1]+(18758.11024597094*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0]+18758.11024597094*rdxLxSq[0])*rdxLy[1]+233.8268590217983*rdxUxCu[0]+15900.22641348229*rdxLx[0]*rdxUxSq[0]+15900.22641348229*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-181.8653347947321*rdxUx[0]*rdxUyCu[1])+((-12003.11209645232*rdxUx[0]*rdxLy[1])+9145.22826396367*rdxUxSq[0]-17874.76433411081*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(12003.11209645232*rdxUx[0]*rdxLySq[1]+9327.093598758403*rdxUxCu[0]+3455.441361099909*rdxLx[0]*rdxUxSq[0]-17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+181.8653347947321*rdxUx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxUxSq[0])*rdxLySq[1]+((-9327.093598758403*rdxUxCu[0])-3455.441361099909*rdxLx[0]*rdxUxSq[0]+17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((-779.4228634059946*rdxLy[1]*rdxUyCu[1])+(((-19303.70625035514*rdxUx[0])-19303.70625035514*rdxLx[0])*rdxLy[1]-52637.02404201817*rdxLySq[1])*rdxUySq[1]+((-28630.79984911354*rdxLyCu[1])+((-68719.1157902952*rdxUx[0])-68719.1157902952*rdxLx[0])*rdxLySq[1]+((-18758.11024597094*rdxUxSq[0])-40893.71956670118*rdxLx[0]*rdxUx[0]-18758.11024597094*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-415.6921938165305*rdxLyR4[1]+((-10574.170180208*rdxUx[0])-10574.170180208*rdxLx[0])*rdxLyCu[1]+((-10392.30484541326*rdxUxSq[0])-47804.60228890101*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-233.8268590217983*rdxUxCu[0])-15900.22641348229*rdxLx[0]*rdxUxSq[0]-15900.22641348229*rdxLxSq[0]*rdxUx[0]-233.8268590217983*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-181.8653347947321*rdxLx[0]*rdxUyCu[1])+((-12003.11209645232*rdxLx[0]*rdxLy[1])-17874.76433411081*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(12003.11209645232*rdxLx[0]*rdxLySq[1]-17692.89899931608*rdxLx[0]*rdxUxSq[0]+3455.441361099909*rdxLxSq[0]*rdxUx[0]+9327.093598758403*rdxLxCu[0])*rdxUy[1]+181.8653347947321*rdxLx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(17692.89899931608*rdxLx[0]*rdxUxSq[0]-3455.441361099909*rdxLxSq[0]*rdxUx[0]-9327.093598758403*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((24375.0*rdxLx[0]-24375.0*rdxUx[0])*rdxUyCu[1]+((32565.0*rdxLx[0]-32565.0*rdxUx[0])*rdxLy[1]-24960.0*rdxUxSq[0]+24960.0*rdxLxSq[0])*rdxUySq[1]+((46215.0*rdxLx[0]-46215.0*rdxUx[0])*rdxLySq[1]+(46800.0*rdxLxSq[0]-46800.0*rdxUxSq[0])*rdxLy[1]-585.0*rdxUxCu[0]-38610.0*rdxLx[0]*rdxUxSq[0]+38610.0*rdxLxSq[0]*rdxUx[0]+585.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]-8640.0*rdxUxSq[0]+19440.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])-8865.0*rdxUxCu[0]-4275.0*rdxLx[0]*rdxUxSq[0]+19215.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+(8640.0*rdxUxSq[0]-19440.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8865.0*rdxUxCu[0]+4275.0*rdxLx[0]*rdxUxSq[0]-19215.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+((46215.0*rdxUx[0]-46215.0*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((32565.0*rdxUx[0]-32565.0*rdxLx[0])*rdxLySq[1]+(46800.0*rdxUxSq[0]-46800.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(24375.0*rdxUx[0]-24375.0*rdxLx[0])*rdxLyCu[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLySq[1]+(585.0*rdxUxCu[0]+38610.0*rdxLx[0]*rdxUxSq[0]-38610.0*rdxLxSq[0]*rdxUx[0]-585.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-19440.0*rdxLx[0]*rdxUx[0]+8640.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-19215.0*rdxLx[0]*rdxUxSq[0]+4275.0*rdxLxSq[0]*rdxUx[0]+8865.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(19440.0*rdxLx[0]*rdxUx[0]-8640.0*rdxLxSq[0])*rdxLySq[1]+(19215.0*rdxLx[0]*rdxUxSq[0]-4275.0*rdxLxSq[0]*rdxUx[0]-8865.0*rdxLxCu[0])*rdxLy[1])*phiLx[2]-396.0*phiUy[1]*rdxUyR4[1]+(((-27378.0*phiUy[1])-846.0*phiLy[1])*rdxLy[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiUy[1]+483.0*rdxUx[0]*phiUx[1]+483.0*rdxLx[0]*phiLx[1]+(23169.64365284887*phiUy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiUy[0])*rdxLx[0])*rdxUyCu[1]+(((-57078.0*phiUy[1])-57078.0*phiLy[1])*rdxLySq[1]+(((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(33842.54072908828*phiUy[0]-40633.91194556586*phiUx[0]+50189.63625092335*phiLy[0])*rdxUx[0]+((-33842.54072908828*phiUy[0])-50189.63625092335*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLy[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiUy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(23777.59348630554*phiUy[0]+21740.70173660454*phiUx[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiUy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-846.0*phiUy[1])-27378.0*phiLy[1])*rdxLyCu[1]+(((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(50189.63625092335*phiUy[0]-40633.91194556586*phiUx[0]+33842.54072908828*phiLy[0])*rdxUx[0]+((-50189.63625092335*phiUy[0])-33842.54072908828*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLySq[1]+(((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiUy[1]+(22980.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0])*phiUx[1]+((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiLy[1]+(88110.0*rdxLx[0]*rdxUx[0]+22980.0*rdxLxSq[0])*phiLx[1]+(50797.58608438003*phiUy[0]-34876.57506120691*phiUx[0]+50797.58608438003*phiLy[0])*rdxUxSq[0]+(102561.6565193835*phiLx[0]-102561.6565193835*phiUx[0])*rdxLx[0]*rdxUx[0]+((-50797.58608438003*phiUy[0])-50797.58608438003*phiLy[0]+34876.57506120691*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiUy[1]+((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(607.9498334566757*phiUy[0]+22712.38223965068*phiUx[0])*rdxUxCu[0]+(40124.68900814059*phiUy[0]-44349.16092780109*phiUx[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40124.68900814059*phiUy[0])-51862.79733103487*phiUx[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiUy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxUy[1]-396.0*phiLy[1]*rdxLyR4[1]+(483.0*rdxUx[0]*phiUx[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiLy[1]+483.0*rdxLx[0]*phiLx[1]+(23169.64365284887*phiLy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiLy[0])*rdxLx[0])*rdxLyCu[1]+((47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiLy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(21740.70173660454*phiUx[0]+23777.59348630554*phiLy[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiLy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiLy[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(22712.38223965068*phiUx[0]+607.9498334566757*phiLy[0])*rdxUxCu[0]+((-44349.16092780109*phiUx[0])+40124.68900814059*phiLy[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-51862.79733103487*phiUx[0])-40124.68900814059*phiLy[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiLy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-396.0*rdxUxR4[0])-25758.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]+774.0*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(774.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]-25758.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[1]+374.1229744348773*phiUx[0]*rdxUxR4[0]+(24224.46259465831*phiUx[0]+841.7766924784738*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(56024.91542162288*phiLx[0]-56024.91542162288*phiUx[0])*rdxLxSq[0]*rdxUxSq[0]+((-841.7766924784738*phiUx[0])-24224.46259465831*phiLx[0])*rdxLxCu[0]*rdxUx[0]-374.1229744348773*phiLx[0]*rdxLxR4[0]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[2] = -(1.0*((((9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxUySq[1]+((7898.15168251408*rdxLx[0]-7898.15168251408*rdxUx[0])*rdxLy[1]+9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxUy[1]+(9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxLySq[1]+(9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxLy[1]+415.6921938165305*rdxUxCu[0]+27435.68479189101*rdxLx[0]*rdxUxSq[0]-27435.68479189101*rdxLxSq[0]*rdxUx[0]-415.6921938165305*rdxLxCu[0])*rho[3]+((-432.0*rdxUyCu[1])+((-29376.0*rdxLy[1])-26928.0*rdxUx[0]-26928.0*rdxLx[0])*rdxUySq[1]+((-29376.0*rdxLySq[1])+((-81936.0*rdxUx[0])-81936.0*rdxLx[0])*rdxLy[1]-27600.0*rdxUxSq[0]-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxUy[1]-432.0*rdxLyCu[1]+((-26928.0*rdxUx[0])-26928.0*rdxLx[0])*rdxLySq[1]+((-27600.0*rdxUxSq[0])-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxLy[1]-1104.0*rdxUxCu[0]-75072.0*rdxLx[0]*rdxUxSq[0]-75072.0*rdxLxSq[0]*rdxUx[0]-1104.0*rdxLxCu[0])*rho[2]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[1]+1080.799703922979*rho[0]*rdxUyCu[1]+(71332.78045891662*rho[0]*rdxLy[1]+(66289.04850727606*rdxUx[0]+66289.04850727606*rdxLx[0])*rho[0])*rdxUySq[1]+((65208.24880335309*rdxUxSq[0]+60164.51685171252*rdxLx[0]*rdxUx[0]+65208.24880335309*rdxLxSq[0])*rho[0]-71332.78045891662*rho[0]*rdxLySq[1])*rdxUy[1]-1080.799703922979*rho[0]*rdxLyCu[1]+((-66289.04850727606*rdxUx[0])-66289.04850727606*rdxLx[0])*rho[0]*rdxLySq[1]+((-65208.24880335309*rdxUxSq[0])-60164.51685171252*rdxLx[0]*rdxUx[0]-65208.24880335309*rdxLxSq[0])*rho[0]*rdxLy[1])*volFac+((9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxUyCu[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLy[1]+9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1]-181.8653347947321*rdxUxCu[0]-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(233.8268590217983*rdxUx[0]*rdxUyCu[1]+(15900.22641348229*rdxUx[0]*rdxLy[1]+10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15900.22641348229*rdxUx[0]*rdxLySq[1]+(47804.60228890101*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0])*rdxLy[1]+10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+415.6921938165305*rdxUxR4[0]+28630.79984911354*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+((-181.8653347947321*rdxUxCu[0])-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15900.22641348229*rdxLx[0]*rdxLy[1])-18758.11024597094*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxUySq[1]+((-15900.22641348229*rdxLx[0]*rdxLySq[1])+((-40893.71956670118*rdxLx[0]*rdxUx[0])-47804.60228890101*rdxLxSq[0])*rdxLy[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxUy[1]-233.8268590217983*rdxLx[0]*rdxLyCu[1]+((-18758.11024597094*rdxLx[0]*rdxUx[0])-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-19303.70625035514*rdxLx[0]*rdxUxSq[0])-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxLy[1]-779.4228634059946*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-28630.79984911354*rdxLxCu[0]*rdxUx[0]-415.6921938165305*rdxLxR4[0])*phiLx[3]+((-396.0*rdxUyR4[1])+((-25758.0*rdxLy[1])-24099.0*rdxUx[0]-24099.0*rdxLx[0])*rdxUyCu[1]+(51462.0*rdxLySq[1]+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLy[1]-23220.0*rdxUxSq[0]+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxUySq[1]+(774.0*rdxLyCu[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1]+483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-243.0*rdxUx[0]*rdxUyCu[1])+((-16524.0*rdxUx[0]*rdxLy[1])-9972.0*rdxUxSq[0]-20322.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-16524.0*rdxUx[0]*rdxLySq[1])+((-50364.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0])*rdxLy[1]-10125.0*rdxUxCu[0]-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-243.0*rdxUx[0]*rdxLyCu[1]+((-9972.0*rdxUxSq[0])-20322.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-10125.0*rdxUxCu[0])-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-396.0*rdxUxR4[0]-27378.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-846.0*rdxLxCu[0]*rdxUx[0])*phiUx[2]+(774.0*rdxLy[1]*rdxUyCu[1]+(51462.0*rdxLySq[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLy[1])*rdxUySq[1]+((-25758.0*rdxLyCu[1])+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-396.0*rdxLyR4[1]+((-24099.0*rdxUx[0])-24099.0*rdxLx[0])*rdxLyCu[1]+((-23220.0*rdxUxSq[0])+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxLySq[1]+(483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-243.0*rdxLx[0]*rdxUyCu[1])+((-16524.0*rdxLx[0]*rdxLy[1])-20322.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*rdxUySq[1]+((-16524.0*rdxLx[0]*rdxLySq[1])+((-41814.0*rdxLx[0]*rdxUx[0])-50364.0*rdxLxSq[0])*rdxLy[1]-20925.0*rdxLx[0]*rdxUxSq[0]-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxUy[1]-243.0*rdxLx[0]*rdxLyCu[1]+((-20322.0*rdxLx[0]*rdxUx[0])-9972.0*rdxLxSq[0])*rdxLySq[1]+((-20925.0*rdxLx[0]*rdxUxSq[0])-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxLy[1]-846.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-27378.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[2]+374.1229744348773*phiUy[0]*rdxUyR4[1]+((24224.46259465831*phiUy[0]+841.7766924784738*phiLy[0])*rdxLy[1]+(8865.0*rdxLx[0]-8865.0*rdxUx[0])*phiUy[1]-585.0*rdxUx[0]*phiUx[1]+585.0*rdxLx[0]*phiLx[1]+(22712.38223965068*phiUy[0]+607.9498334566757*phiUx[0])*rdxUx[0]+(22712.38223965068*phiUy[0]+607.9498334566757*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((56024.91542162288*phiLy[0]-56024.91542162288*phiUy[0])*rdxLySq[1]+((4275.0*rdxLx[0]-4275.0*rdxUx[0])*phiUy[1]-38610.0*rdxUx[0]*phiUx[1]+(19215.0*rdxLx[0]-19215.0*rdxUx[0])*phiLy[1]+38610.0*rdxLx[0]*phiLx[1]+((-44349.16092780109*phiUy[0])+40124.68900814059*phiUx[0]+51862.79733103487*phiLy[0])*rdxUx[0]+((-44349.16092780109*phiUy[0])+51862.79733103487*phiLy[0]+40124.68900814059*phiLx[0])*rdxLx[0])*rdxLy[1]+(8640.0*rdxLxSq[0]-8640.0*rdxUxSq[0])*phiUy[1]+((-24960.0*rdxUxSq[0])-46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(46800.0*rdxLx[0]*rdxUx[0]+24960.0*rdxLxSq[0])*phiLx[1]+(21740.70173660454*phiUy[0]+23777.59348630554*phiUx[0])*rdxUxSq[0]+((-34876.57506120691*phiUy[0])+50797.58608438003*phiUx[0]+50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+(21740.70173660454*phiUy[0]+23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-841.7766924784738*phiUy[0])-24224.46259465831*phiLy[0])*rdxLyCu[1]+((19215.0*rdxUx[0]-19215.0*rdxLx[0])*phiUy[1]+38610.0*rdxUx[0]*phiUx[1]+(4275.0*rdxUx[0]-4275.0*rdxLx[0])*phiLy[1]-38610.0*rdxLx[0]*phiLx[1]+((-51862.79733103487*phiUy[0])-40124.68900814059*phiUx[0]+44349.16092780109*phiLy[0])*rdxUx[0]+((-51862.79733103487*phiUy[0])+44349.16092780109*phiLy[0]-40124.68900814059*phiLx[0])*rdxLx[0])*rdxLySq[1]+((19440.0*rdxUxSq[0]-19440.0*rdxLxSq[0])*phiUy[1]+(19440.0*rdxLxSq[0]-19440.0*rdxUxSq[0])*phiLy[1]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxUxSq[0]+(102561.6565193835*phiLy[0]-102561.6565193835*phiUy[0])*rdxLx[0]*rdxUx[0]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxLxSq[0])*rdxLy[1]+(225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*phiUy[1]+((-24375.0*rdxUxCu[0])-32565.0*rdxLx[0]*rdxUxSq[0]-46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(46215.0*rdxLx[0]*rdxUxSq[0]+32565.0*rdxLxSq[0]*rdxUx[0]+24375.0*rdxLxCu[0])*phiLx[1]+(23169.64365284887*phiUx[0]-597.5575286112626*phiUy[0])*rdxUxCu[0]+((-40633.91194556586*phiUy[0])+33842.54072908828*phiUx[0]+50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40633.91194556586*phiUy[0])+50189.63625092335*phiUx[0]+33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(23169.64365284887*phiLx[0]-597.5575286112626*phiUy[0])*rdxLxCu[0])*rdxUy[1]-374.1229744348773*phiLy[0]*rdxLyR4[1]+(585.0*rdxUx[0]*phiUx[1]+(8865.0*rdxUx[0]-8865.0*rdxLx[0])*phiLy[1]-585.0*rdxLx[0]*phiLx[1]+((-607.9498334566757*phiUx[0])-22712.38223965068*phiLy[0])*rdxUx[0]+((-22712.38223965068*phiLy[0])-607.9498334566757*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((24960.0*rdxUxSq[0]+46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(8640.0*rdxUxSq[0]-8640.0*rdxLxSq[0])*phiLy[1]+((-46800.0*rdxLx[0]*rdxUx[0])-24960.0*rdxLxSq[0])*phiLx[1]+((-23777.59348630554*phiUx[0])-21740.70173660454*phiLy[0])*rdxUxSq[0]+((-50797.58608438003*phiUx[0])+34876.57506120691*phiLy[0]-50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+((-21740.70173660454*phiLy[0])-23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((24375.0*rdxUxCu[0]+32565.0*rdxLx[0]*rdxUxSq[0]+46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*phiLy[1]+((-46215.0*rdxLx[0]*rdxUxSq[0])-32565.0*rdxLxSq[0]*rdxUx[0]-24375.0*rdxLxCu[0])*phiLx[1]+(597.5575286112626*phiLy[0]-23169.64365284887*phiUx[0])*rdxUxCu[0]+((-33842.54072908828*phiUx[0])+40633.91194556586*phiLy[0]-50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-50189.63625092335*phiUx[0])+40633.91194556586*phiLy[0]-33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(597.5575286112626*phiLy[0]-23169.64365284887*phiLx[0])*rdxLxCu[0])*rdxLy[1]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[3] = (((144.0*rdxUyCu[1]+(9792.0*rdxLy[1]+3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxUySq[1]+(9792.0*rdxLySq[1]+(31568.0*rdxUx[0]+31568.0*rdxLx[0])*rdxLy[1]+3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxUy[1]+144.0*rdxLyCu[1]+(3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxLySq[1]+(3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxLy[1]+144.0*rdxUxCu[0]+9792.0*rdxLx[0]*rdxUxSq[0]+9792.0*rdxLxSq[0]*rdxUx[0]+144.0*rdxLxCu[0])*rho[3]+((8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxUySq[1]+((6845.064791512203*rdxUx[0]-6845.064791512203*rdxLx[0])*rdxLy[1]-8646.397631383834*rdxUxSq[0]+8646.397631383834*rdxLxSq[0])*rdxUy[1]+(8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxLySq[1]+(8646.397631383834*rdxLxSq[0]-8646.397631383834*rdxUxSq[0])*rdxLy[1]-360.2665679743264*rdxUxCu[0]-23777.59348630554*rdxLx[0]*rdxUxSq[0]+23777.59348630554*rdxLxSq[0]*rdxUx[0]+360.2665679743264*rdxLxCu[0])*rho[2]+((-360.2665679743264*rdxUyCu[1])+((-23777.59348630554*rdxLy[1])-8646.397631383834*rdxUx[0]-8646.397631383834*rdxLx[0])*rdxUySq[1]+(23777.59348630554*rdxLySq[1]-8286.131063409508*rdxUxSq[0]+6845.064791512203*rdxLx[0]*rdxUx[0]-8286.131063409508*rdxLxSq[0])*rdxUy[1]+360.2665679743264*rdxLyCu[1]+(8646.397631383834*rdxUx[0]+8646.397631383834*rdxLx[0])*rdxLySq[1]+(8286.131063409508*rdxUxSq[0]-6845.064791512203*rdxLx[0]*rdxUx[0]+8286.131063409508*rdxLxSq[0])*rdxLy[1])*rho[1]+(21632.0*rdxUx[0]-21632.0*rdxLx[0])*rho[0]*rdxUySq[1]+(21632.0*rdxUxSq[0]-21632.0*rdxLxSq[0])*rho[0]*rdxUy[1]+(21632.0*rdxLx[0]-21632.0*rdxUx[0])*rho[0]*rdxLySq[1]+(21632.0*rdxLxSq[0]-21632.0*rdxUxSq[0])*rho[0]*rdxLy[1])*volFac+(132.0*rdxUyR4[1]+(8586.0*rdxLy[1]+3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxUyCu[1]+((-17154.0*rdxLySq[1])+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLy[1]+2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-258.0*rdxLyCu[1])+((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1]-63.0*rdxUxCu[0]-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-63.0*rdxUx[0]*rdxUyCu[1])+((-4284.0*rdxUx[0]*rdxLy[1])+2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-4284.0*rdxUx[0]*rdxLySq[1])+((-17516.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0])*rdxLy[1]+3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-63.0*rdxUx[0]*rdxLyCu[1]+(2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+132.0*rdxUxR4[0]+8586.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]-258.0*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((-258.0*rdxLy[1]*rdxUyCu[1])+(((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLy[1]-17154.0*rdxLySq[1])*rdxUySq[1]+(8586.0*rdxLyCu[1]+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+132.0*rdxLyR4[1]+(3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxLyCu[1]+(2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxLySq[1]+((-63.0*rdxUxCu[0])-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-63.0*rdxLx[0]*rdxUyCu[1])+((-4284.0*rdxLx[0]*rdxLy[1])-6158.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-4284.0*rdxLx[0]*rdxLySq[1])+((-10106.0*rdxLx[0]*rdxUx[0])-17516.0*rdxLxSq[0])*rdxLy[1]-6353.0*rdxLx[0]*rdxUxSq[0]-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxUy[1]-63.0*rdxLx[0]*rdxLyCu[1]+(2812.0*rdxLxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6353.0*rdxLx[0]*rdxUxSq[0])-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxLy[1]-258.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]+8586.0*rdxLxCu[0]*rdxUx[0]+132.0*rdxLxR4[0])*phiLx[3]+((8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxUyCu[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLy[1]-7925.864495435182*rdxUxSq[0]+7925.864495435182*rdxLxSq[0])*rdxUySq[1]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1]+157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(77.94228634059945*rdxUx[0]*rdxUyCu[1]+(5300.075471160763*rdxUx[0]*rdxLy[1]-2591.14800812304*rdxUxSq[0]+6730.749438212657*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(5300.075471160763*rdxUx[0]*rdxLySq[1]+(20937.03016189259*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0])*rdxLy[1]-2793.797952608599*rdxUxCu[0]+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+77.94228634059945*rdxUx[0]*rdxLyCu[1]+(6730.749438212657*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxUxSq[0])*rdxLySq[1]+((-2793.797952608599*rdxUxCu[0])+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-124.7076581449591*rdxUxR4[0]-8074.820864886104*rdxLx[0]*rdxUxCu[0]+18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+280.592230826158*rdxLxCu[0]*rdxUx[0])*phiUx[2]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxLyCu[1]+(7925.864495435182*rdxLxSq[0]-7925.864495435182*rdxUxSq[0])*rdxLySq[1]+(157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-77.94228634059945*rdxLx[0]*rdxUyCu[1])+((-5300.075471160763*rdxLx[0]*rdxLy[1])-6730.749438212657*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*rdxUySq[1]+((-5300.075471160763*rdxLx[0]*rdxLySq[1])+((-13236.33227144136*rdxLx[0]*rdxUx[0])-20937.03016189259*rdxLxSq[0])*rdxLy[1]-6933.399382698215*rdxLx[0]*rdxUxSq[0]-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxUy[1]-77.94228634059945*rdxLx[0]*rdxLyCu[1]+(2591.14800812304*rdxLxSq[0]-6730.749438212657*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6933.399382698215*rdxLx[0]*rdxUxSq[0])-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxLy[1]-280.592230826158*rdxLx[0]*rdxUxCu[0]-18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+8074.820864886104*rdxLxCu[0]*rdxUx[0]+124.7076581449591*rdxLxR4[0])*phiLx[2]-124.7076581449591*phiUy[1]*rdxUyR4[1]+(((-8074.820864886104*phiUy[1])-280.592230826158*phiLy[1])*rdxLy[1]+((-2793.797952608599*rdxUx[0])-2793.797952608599*rdxLx[0])*phiUy[1]+157.6166234887678*rdxUx[0]*phiUx[1]+157.6166234887678*rdxLx[0]*phiLx[1]+(7683.0*phiUy[0]-195.0*phiUx[0])*rdxUx[0]+(195.0*phiLx[0]-7683.0*phiUy[0])*rdxLx[0])*rdxUyCu[1]+((18674.97180720763*phiUy[1]-18674.97180720763*phiLy[1])*rdxLySq[1]+((17086.68121666697*rdxUx[0]+17086.68121666697*rdxLx[0])*phiUy[1]+10402.69715025867*rdxUx[0]*phiUx[1]+((-6933.399382698215*rdxUx[0])-6933.399382698215*rdxLx[0])*phiLy[1]+10402.69715025867*rdxLx[0]*phiLx[1]+(3705.0*phiUy[0]-12870.0*phiUx[0]+16653.0*phiLy[0])*rdxUx[0]+((-3705.0*phiUy[0])-16653.0*phiLy[0]+12870.0*phiLx[0])*rdxLx[0])*rdxLy[1]+((-2591.14800812304*rdxUxSq[0])+20937.03016189259*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxLxSq[0])*phiUy[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxUxSq[0])*phiUx[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxLxSq[0])*phiLx[1]+(7488.0*phiUy[0]+7488.0*phiUx[0])*rdxUxSq[0]+(16848.0*phiLx[0]-16848.0*phiUx[0])*rdxLx[0]*rdxUx[0]+((-7488.0*phiUy[0])-7488.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((280.592230826158*phiUy[1]+8074.820864886104*phiLy[1])*rdxLyCu[1]+((6933.399382698215*rdxUx[0]+6933.399382698215*rdxLx[0])*phiUy[1]-10402.69715025867*rdxUx[0]*phiUx[1]+((-17086.68121666697*rdxUx[0])-17086.68121666697*rdxLx[0])*phiLy[1]-10402.69715025867*rdxLx[0]*phiLx[1]+((-16653.0*phiUy[0])+12870.0*phiUx[0]-3705.0*phiLy[0])*rdxUx[0]+(16653.0*phiUy[0]+3705.0*phiLy[0]-12870.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((6730.749438212657*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0]+6730.749438212657*rdxLxSq[0])*phiUy[1]+((-6730.749438212657*rdxUxSq[0])-13236.33227144136*rdxLx[0]*rdxUx[0]-6730.749438212657*rdxLxSq[0])*phiLy[1]+(16848.0*phiLy[0]-16848.0*phiUy[0])*rdxUxSq[0]+(16848.0*phiUy[0]-16848.0*phiLy[0])*rdxLxSq[0])*rdxLy[1]+(77.94228634059945*rdxUxCu[0]+5300.075471160763*rdxLx[0]*rdxUxSq[0]+5300.075471160763*rdxLxSq[0]*rdxUx[0]+77.94228634059945*rdxLxCu[0])*phiUy[1]+((-8083.48111892395*rdxUxCu[0])-2994.715846286589*rdxLx[0]*rdxUxSq[0]+15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(15333.84579940727*rdxLx[0]*rdxUxSq[0]-2994.715846286589*rdxLxSq[0]*rdxUx[0]-8083.48111892395*rdxLxCu[0])*phiLx[1]+(7683.0*phiUx[0]-195.0*phiUy[0])*rdxUxCu[0]+((-12870.0*phiUy[0])+3705.0*phiUx[0]+16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(12870.0*phiUy[0]-16653.0*phiUx[0]-3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(195.0*phiUy[0]-7683.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+124.7076581449591*phiLy[1]*rdxLyR4[1]+((-157.6166234887678*rdxUx[0]*phiUx[1])+(2793.797952608599*rdxUx[0]+2793.797952608599*rdxLx[0])*phiLy[1]-157.6166234887678*rdxLx[0]*phiLx[1]+(195.0*phiUx[0]-7683.0*phiLy[0])*rdxUx[0]+(7683.0*phiLy[0]-195.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((7925.864495435182*rdxUxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiUx[1]+(2591.14800812304*rdxUxSq[0]-20937.03016189259*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*phiLy[1]+(7925.864495435182*rdxLxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiLx[1]+((-7488.0*phiUx[0])-7488.0*phiLy[0])*rdxUxSq[0]+(16848.0*phiUx[0]-16848.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(7488.0*phiLy[0]+7488.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((8083.48111892395*rdxUxCu[0]+2994.715846286589*rdxLx[0]*rdxUxSq[0]-15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-77.94228634059945*rdxUxCu[0])-5300.075471160763*rdxLx[0]*rdxUxSq[0]-5300.075471160763*rdxLxSq[0]*rdxUx[0]-77.94228634059945*rdxLxCu[0])*phiLy[1]+((-15333.84579940727*rdxLx[0]*rdxUxSq[0])+2994.715846286589*rdxLxSq[0]*rdxUx[0]+8083.48111892395*rdxLxCu[0])*phiLx[1]+(195.0*phiLy[0]-7683.0*phiUx[0])*rdxUxCu[0]+((-3705.0*phiUx[0])+12870.0*phiLy[0]-16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(16653.0*phiUx[0]-12870.0*phiLy[0]+3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(7683.0*phiLx[0]-195.0*phiLy[0])*rdxLxCu[0])*rdxLy[1])/(12.0*rdxUyR4[1]+(1608.0*rdxLy[1]+1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxUyCu[1]+(53892.0*rdxLySq[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLy[1]+2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxUySq[1]+(1608.0*rdxLyCu[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLySq[1]+(69048.0*rdxUxSq[0]+166696.0*rdxLx[0]*rdxUx[0]+69048.0*rdxLxSq[0])*rdxLy[1]+1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxUy[1]+12.0*rdxLyR4[1]+(1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxLyCu[1]+(2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxLySq[1]+(1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxLy[1]+12.0*rdxUxR4[0]+1608.0*rdxLx[0]*rdxUxCu[0]+53892.0*rdxLxSq[0]*rdxUxSq[0]+1608.0*rdxLxCu[0]*rdxUx[0]+12.0*rdxLxR4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[0]*rho[1]+144.0*rho[0]*rdxCp2[1]+1600.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(81.0*phiUy[0]+81.0*phiLy[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+93.53074360871933*rdxCp2[0]*phiUx[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+(900.0*phiUy[0]-54.0*phiUx[0]+900.0*phiLy[0]+864.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiUx[1]+(1020.0*phiUx[0]+3120.0*bcVals[0])*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]+277.1281292110203*rdxCp2[0]*rho[0])*volFac+((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]-360.0*rdxCp2[0]*phiUx[1]+189.0*rdxCp2[0]*phiLy[1]+(155.8845726811989*phiUy[0]+311.7691453623978*phiUx[0]+155.8845726811989*phiLy[0]-1247.076581449591*bcVals[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiUx[1]+(623.5382907247956*phiUx[0]-1247.076581449591*bcVals[0])*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = ((748.2459488697547*rdxCp2[0]*rho[3]+(368.0*rdxCp2[1]+1600.0*rdxCp2[0])*rho[2])*volFac-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(1020.0*rdxCp2Sq[0]-138.0*rdxCp2[0]*rdxCp2[1])*phiUx[2]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(199.1858428704209*phiUy[0]-199.1858428704209*phiLy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiUy[0]-866.0254037844386*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[0]*rho[2])*volFac+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiUx[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+623.5382907247956*rdxCp2Sq[0])*phiUx[2]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiUy[0]-150.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[1]-576.0*rho[0]*rdxCp2[1]-1680.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-324.0*phiUy[0])-324.0*phiLy[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+311.7691453623978*rdxCp2[0]*phiUx[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+((-945.0*phiUy[0])-324.0*phiUx[0]-945.0*phiLy[0]+576.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiUx[1]+(2080.0*bcVals[0]-720.0*phiUx[0])*rdxCp2Sq[0]))/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]-138.5640646055102*rdxCp2[0]*rho[0])*volFac+((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]+75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiUx[1]+54.0*rdxCp2[0]*phiLy[1]+((-77.94228634059945*phiUy[0])+155.8845726811989*phiUx[0]-77.94228634059945*phiLy[0]+277.1281292110203*bcVals[0])*rdxCp2[0])*rdxCp2[1]+277.1281292110203*bcVals[0]*rdxCp2Sq[0])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[3]+((-1472.0*rdxCp2[1])-1680.0*rdxCp2[0])*rho[2])*volFac-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-828.0*rdxCp2[0]*rdxCp2[1])-720.0*rdxCp2Sq[0])*phiUx[2]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(796.7433714816835*phiLy[0]-796.7433714816835*phiUy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiLy[0]-909.3266739736605*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[0]*rho[2])*volFac+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]+1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiUx[2]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiLy[0]-225.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[1]-144.0*rho[0]*rdxCp2[1]-1600.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-81.0*phiUy[0])-81.0*phiLy[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+93.53074360871933*rdxCp2[0]*phiLx[1]-864.0*rdxCp2[0]*bcVals[1]+((-900.0*phiUy[0])-900.0*phiLy[0]+54.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiLx[1]-3120.0*rdxCp2Sq[0]*bcVals[1]-1020.0*phiLx[0]*rdxCp2Sq[0]))/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]-277.1281292110203*rdxCp2[0]*rho[0])*volFac+((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]+189.0*rdxCp2[0]*phiLy[1]-360.0*rdxCp2[0]*phiLx[1]+1247.076581449591*rdxCp2[0]*bcVals[1]+((-155.8845726811989*phiUy[0])-155.8845726811989*phiLy[0]-311.7691453623978*phiLx[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiLx[1]+1247.076581449591*rdxCp2Sq[0]*bcVals[1]-623.5382907247956*phiLx[0]*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[3]+((-368.0*rdxCp2[1])-1600.0*rdxCp2[0])*rho[2])*volFac-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiLx[3]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(138.0*rdxCp2[0]*rdxCp2[1]-1020.0*rdxCp2Sq[0])*phiLx[2]+(199.1858428704209*phiLy[0]-199.1858428704209*phiUy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiLy[0]-866.0254037844386*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[0]*rho[2])*volFac+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiLx[3]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-796.7433714816835*rdxCp2[0]*rdxCp2[1])-623.5382907247956*rdxCp2Sq[0])*phiLx[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiLy[0]-150.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[0]*rho[1]+576.0*rho[0]*rdxCp2[1]+1680.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(324.0*phiUy[0]+324.0*phiLy[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+311.7691453623978*rdxCp2[0]*phiLx[1]+576.0*rdxCp2[0]*bcVals[1]+(945.0*phiUy[0]+945.0*phiLy[0]+324.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiLx[1]+2080.0*rdxCp2Sq[0]*bcVals[1]+720.0*phiLx[0]*rdxCp2Sq[0])/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]+138.5640646055102*rdxCp2[0]*rho[0])*volFac+((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]-75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]+54.0*rdxCp2[0]*phiLy[1]-150.0*rdxCp2[0]*phiLx[1]+277.1281292110203*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiUy[0]+77.94228634059945*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1]+277.1281292110203*rdxCp2Sq[0]*bcVals[1])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = ((277.1281292110203*rdxCp2[0]*rho[3]+(1472.0*rdxCp2[1]+1680.0*rdxCp2[0])*rho[2])*volFac-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiLx[3]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(828.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiLx[2]+(796.7433714816835*phiUy[0]-796.7433714816835*phiLy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiUy[0]-909.3266739736605*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[0]*rho[2])*volFac+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]-1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiUy[0]-225.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[1]*rho[2]+1600.0*rho[0]*rdxCp2[1]+144.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(3120.0*rdxCp2Sq[1]+864.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+1020.0*phiUy[0]*rdxCp2Sq[1]+((-866.0254037844386*rdxCp2[0]*phiUx[1])+866.0254037844386*rdxCp2[0]*phiLx[1]+((-54.0*phiUy[0])+900.0*phiUx[0]+900.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-77.94228634059945*rdxCp2Sq[0]*phiUx[1]+77.94228634059945*rdxCp2Sq[0]*phiLx[1]+(81.0*phiUx[0]+81.0*phiLx[0])*rdxCp2Sq[0])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = ((748.2459488697547*rdxCp2[1]*rho[3]+(1600.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[1])*volFac+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+1020.0*phiUy[1]*rdxCp2Sq[1]+((-138.0*rdxCp2[0]*phiUy[1])-700.0*rdxCp2[0]*phiUx[1]-700.0*rdxCp2[0]*phiLx[1]+(866.0254037844386*phiUx[0]-866.0254037844386*phiLx[0])*rdxCp2[0])*rdxCp2[1]-161.0*rdxCp2Sq[0]*phiUx[1]-161.0*rdxCp2Sq[0]*phiLx[1]+(199.1858428704209*phiUx[0]-199.1858428704209*phiLx[0])*rdxCp2Sq[0])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]+277.1281292110203*rho[0]*rdxCp2[1])*volFac+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]+((-1247.076581449591*rdxCp2Sq[1])-1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[2]+623.5382907247956*phiUy[0]*rdxCp2Sq[1]+((-150.0*rdxCp2[0]*phiUx[1])+150.0*rdxCp2[0]*phiLx[1]+(311.7691453623978*phiUy[0]+155.8845726811989*phiUx[0]+155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[1]*rho[1])*volFac+((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]+623.5382907247956*phiUy[1]*rdxCp2Sq[1]+(796.7433714816835*rdxCp2[0]*phiUy[1]-121.2435565298214*rdxCp2[0]*phiUx[1]-121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiUx[0]-150.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[2]-1680.0*rho[0]*rdxCp2[1]-576.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]-720.0*phiUy[0]*rdxCp2Sq[1]+(909.3266739736605*rdxCp2[0]*phiUx[1]-909.3266739736605*rdxCp2[0]*phiLx[1]+((-324.0*phiUy[0])-945.0*phiUx[0]-945.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+311.7691453623978*rdxCp2Sq[0]*phiUx[1]-311.7691453623978*rdxCp2Sq[0]*phiLx[1]+((-324.0*phiUx[0])-324.0*phiLx[0])*rdxCp2Sq[0]))/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[3]+((-1680.0*rdxCp2[1])-1472.0*rdxCp2[0])*rho[1])*volFac+(692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiUy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-720.0*phiUy[1]*rdxCp2Sq[1]+((-828.0*rdxCp2[0]*phiUy[1])+735.0*rdxCp2[0]*phiUx[1]+735.0*rdxCp2[0]*phiLx[1]+(909.3266739736605*phiLx[0]-909.3266739736605*phiUx[0])*rdxCp2[0])*rdxCp2[1]+644.0*rdxCp2Sq[0]*phiUx[1]+644.0*rdxCp2Sq[0]*phiLx[1]+(796.7433714816835*phiLx[0]-796.7433714816835*phiUx[0])*rdxCp2Sq[0]))/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]-138.5640646055102*rho[0]*rdxCp2[1])*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(75.0*rdxCp2[0]*phiUx[1]-75.0*rdxCp2[0]*phiLx[1]+(155.8845726811989*phiUy[0]-77.94228634059945*phiUx[0]-77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[1]*rho[1])*volFac-1150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+(1195.115057222525*rdxCp2[0]*phiUy[1]+181.8653347947321*rdxCp2[0]*phiUx[1]+181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiLx[0]-225.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[2]-1600.0*rho[0]*rdxCp2[1]-144.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-3120.0*rdxCp2Sq[1])-864.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]-1020.0*phiLy[0]*rdxCp2Sq[1]+(866.0254037844386*rdxCp2[0]*phiUx[1]-866.0254037844386*rdxCp2[0]*phiLx[1]+((-900.0*phiUx[0])+54.0*phiLy[0]-900.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0]*phiUx[1]-77.94228634059945*rdxCp2Sq[0]*phiLx[1]+((-81.0*phiUx[0])-81.0*phiLx[0])*rdxCp2Sq[0]))/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[3]+((-1600.0*rdxCp2[1])-368.0*rdxCp2[0])*rho[1])*volFac-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-1020.0*phiLy[1]*rdxCp2Sq[1]+(700.0*rdxCp2[0]*phiUx[1]+138.0*rdxCp2[0]*phiLy[1]+700.0*rdxCp2[0]*phiLx[1]+(866.0254037844386*phiLx[0]-866.0254037844386*phiUx[0])*rdxCp2[0])*rdxCp2[1]+161.0*rdxCp2Sq[0]*phiUx[1]+161.0*rdxCp2Sq[0]*phiLx[1]+(199.1858428704209*phiLx[0]-199.1858428704209*phiUx[0])*rdxCp2Sq[0]))/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]-277.1281292110203*rho[0]*rdxCp2[1])*volFac+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+(1247.076581449591*rdxCp2Sq[1]+1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]-623.5382907247956*phiLy[0]*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiLx[1]+((-155.8845726811989*phiUx[0])-311.7691453623978*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[1]*rho[1])*volFac+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]-623.5382907247956*phiLy[1]*rdxCp2Sq[1]+(121.2435565298214*rdxCp2[0]*phiUx[1]-796.7433714816835*rdxCp2[0]*phiLy[1]+121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiLx[0]-150.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[1]*rho[2]+1680.0*rho[0]*rdxCp2[1]+576.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+720.0*phiLy[0]*rdxCp2Sq[1]+((-909.3266739736605*rdxCp2[0]*phiUx[1])+909.3266739736605*rdxCp2[0]*phiLx[1]+(945.0*phiUx[0]+324.0*phiLy[0]+945.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-311.7691453623978*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLx[1]+(324.0*phiUx[0]+324.0*phiLx[0])*rdxCp2Sq[0])/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = ((277.1281292110203*rdxCp2[1]*rho[3]+(1680.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[1])*volFac-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3]+(692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiLy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+720.0*phiLy[1]*rdxCp2Sq[1]+((-735.0*rdxCp2[0]*phiUx[1])+828.0*rdxCp2[0]*phiLy[1]-735.0*rdxCp2[0]*phiLx[1]+(909.3266739736605*phiUx[0]-909.3266739736605*phiLx[0])*rdxCp2[0])*rdxCp2[1]-644.0*rdxCp2Sq[0]*phiUx[1]-644.0*rdxCp2Sq[0]*phiLx[1]+(796.7433714816835*phiUx[0]-796.7433714816835*phiLx[0])*rdxCp2Sq[0])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]+138.5640646055102*rho[0]*rdxCp2[1])*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+((-75.0*rdxCp2[0]*phiUx[1])+75.0*rdxCp2[0]*phiLx[1]+(77.94228634059945*phiUx[0]-155.8845726811989*phiLy[0]+77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[1]*rho[1])*volFac+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+((-181.8653347947321*rdxCp2[0]*phiUx[1])-1195.115057222525*rdxCp2[0]*phiLy[1]-181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiUx[0]-225.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+328950.0*phiUy[0]*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiUx[1]+(2715915.0*phiUy[0]-193500.0*phiUx[0]+3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiUx[1]+((-193500.0*phiUy[0])+2715915.0*phiUx[0]+9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiUx[1]+(328950.0*phiUx[0]+1006200.0*bcVals[0])*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+328950.0*phiUy[1]*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiUx[1]+(567939.4598018348*phiUy[0]+1117172.770881926*phiUx[0]-4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiUx[1]+((-33515.18312645776*phiUy[0])+1919718.512568964*phiUx[0]-4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiUx[1]+(201091.0987587466*phiUx[0]-402182.1975174932*bcVals[0])*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiUx[1]+(1919718.512568964*phiUy[0]-33515.18312645776*phiUx[0]+536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiUx[1]+(1117172.770881926*phiUy[0]+567939.4598018348*phiUx[0]+1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiUx[1]+(44400.0*phiUy[0]+64500.0*phiUx[0]-258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiUx[1]+(64500.0*phiUy[0]+44400.0*phiUx[0]+40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*volFac+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-12960.0*phiUy[0]*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiUx[1]+((-456408.0*phiUy[0])+11340.0*phiUx[0]-181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiUx[1]+((-1842246.0*phiUy[0])-532953.0*phiUx[0]-2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiUx[1]+((-928800.0*phiUy[0])-2563956.0*phiUx[0]-8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiUx[1]+((-1052640.0*phiUx[0])-3219840.0*bcVals[0])*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*volFac+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-12960.0*phiUy[1]*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiUx[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiUx[1]+((-301792.5327108011*phiUy[0])-659547.6270141526*phiUx[0]+1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiUx[1]+((-160872.8790069972*phiUy[0])-1745283.675738703*phiUx[0]+3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiUx[1]+(1286983.032055978*bcVals[0]-643491.516027989*phiUx[0])*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiUx[1]+(98207.28078915528*phiUy[0]+14029.6115413079*phiUx[0]-224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiUx[1]+(2061183.762277152*phiUy[0]-814886.6036909672*phiUx[0]-2492594.31717237*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiUx[1]+(6703036.625291553*phiUy[0]-3407636.758811008*phiUx[0]-6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiUx[1]+(332100.0*bcVals[0]-166050.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiUx[1]+(387000.0*phiUy[0]-266400.0*phiUx[0]-241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*volFac+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-1052640.0*phiUy[0]*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiUx[1]+((-2563956.0*phiUy[0])-928800.0*phiUx[0]+1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiUx[1]+((-532953.0*phiUy[0])-1842246.0*phiUx[0]+5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiUx[1]+(11340.0*phiUy[0]-456408.0*phiUx[0]+1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiUx[1]+(37440.0*bcVals[0]-12960.0*phiUx[0])*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+5263200.0*phiUy[1]*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiUx[1]+((-3407636.758811008*phiUy[0])+6703036.625291553*phiUx[0]+1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiUx[1]+((-814886.6036909672*phiUy[0])+2061183.762277152*phiUx[0]+1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiUx[1]+(14029.6115413079*phiUy[0]+98207.28078915528*phiUx[0]+2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+74824.59488697546*bcVals[0]*rdxCp2R4[0])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiUx[1]+((-1745283.675738703*phiUy[0])-160872.8790069972*phiUx[0]+285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiUx[1]+((-659547.6270141526*phiUy[0])-301792.5327108011*phiUx[0]+871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiUx[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiUx[1]+((-266400.0*phiUy[0])+387000.0*phiUx[0]+688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+(580800.0*bcVals[0]-166050.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*volFac+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+1440.0*phiUy[0]*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiUx[1]+(11799.0*phiUy[0]+1890.0*phiUx[0]-3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiUx[1]+(1890.0*phiUy[0]+11799.0*phiUx[0]-33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiUx[1]+(1440.0*phiUx[0]-4160.0*bcVals[0])*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-1440.0*phiUy[1]*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiUx[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiUx[0]-4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiUx[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiUx[0]-14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1662.768775266122*bcVals[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiUx[1]+((-467.6537180435967*phiUy[0])+467.6537180435967*phiUx[0]-831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiUx[1]+((-2727.980021920981*phiUy[0])+2727.980021920981*phiUx[0]-8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiUx[1]+(1125.0*phiUy[0]-1125.0*phiUx[0]-2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiUx[1]+((-1125.0*phiUy[0])+1125.0*phiUx[0]-4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-328950.0*phiLy[0]*rdxCp2R3[1]+((-335151.8312645776*rdxCp2[0]*phiUx[1])-1533436.541464953*rdxCp2[0]*phiLy[1]+(193500.0*phiUx[0]-2715915.0*phiLy[0]-3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2176434.423012786*rdxCp2Sq[0]*phiUx[1]+90490.99444143593*rdxCp2Sq[0]*phiLy[1]+((-2715915.0*phiUx[0])+193500.0*phiLy[0]-9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiUx[1]+((-328950.0*phiUx[0])-1006200.0*bcVals[0])*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]-328950.0*phiLy[1]*rdxCp2R3[1]+(1290000.0*rdxCp2[0]*phiUx[1]-125505.0*rdxCp2[0]*phiLy[1]+((-1117172.770881926*phiUx[0])-567939.4598018348*phiLy[0]+4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2054550.0*rdxCp2Sq[0]*phiUx[1]+40635.0*rdxCp2Sq[0]*phiLy[1]+((-1919718.512568964*phiUx[0])+33515.18312645776*phiLy[0]+4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiUx[1]+(402182.1975174932*bcVals[0]-201091.0987587466*phiUx[0])*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-58050.0*rdxCp2[0]*phiUx[1])-359640.0*rdxCp2[0]*phiLy[1]+(33515.18312645776*phiUx[0]-1919718.512568964*phiLy[0]-536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(493650.0*rdxCp2Sq[0]*phiUx[1]-522450.0*rdxCp2Sq[0]*phiLy[1]+((-567939.4598018348*phiUx[0])-1117172.770881926*phiLy[0]-1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+(74478.18472546172*rdxCp2[0]*phiUx[1]-437394.7904353685*rdxCp2[0]*phiLy[1]+((-64500.0*phiUx[0])-44400.0*phiLy[0]+258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(40789.79651824706*rdxCp2Sq[0]*phiUx[1]-78202.0939617348*rdxCp2Sq[0]*phiLy[1]+((-44400.0*phiUx[0])-64500.0*phiLy[0]-40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*volFac+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+12960.0*phiLy[0]*rdxCp2R4[1]+(19641.45615783106*rdxCp2[0]*phiUx[1]+176773.1054204795*rdxCp2[0]*phiLy[1]+((-11340.0*phiUx[0])+456408.0*phiLy[0]+181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-386469.0325912283*rdxCp2Sq[0]*phiUx[1])+814839.8383191626*rdxCp2Sq[0]*phiLy[1]+(532953.0*phiUx[0]+1842246.0*phiLy[0]+2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2064285.865273508*rdxCp2R3[0]*phiUx[1])+434356.7733188925*rdxCp2R3[0]*phiLy[1]+(2563956.0*phiUx[0]+928800.0*phiLy[0]+8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiUx[1]+(1052640.0*phiUx[0]+3219840.0*bcVals[0])*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*volFac+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+12960.0*phiLy[1]*rdxCp2R4[1]+((-75600.0*rdxCp2[0]*phiUx[1])+157788.0*rdxCp2[0]*phiLy[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-727155.0*rdxCp2Sq[0]*phiUx[1])+465750.0*rdxCp2Sq[0]*phiLy[1]+(659547.6270141526*phiUx[0]+301792.5327108011*phiLy[0]-1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1862820.0*rdxCp2R3[0]*phiUx[1])+195048.0*rdxCp2R3[0]*phiLy[1]+(1745283.675738703*phiUx[0]+160872.8790069972*phiLy[0]-3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiUx[1]+(643491.516027989*phiUx[0]-1286983.032055978*bcVals[0])*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(24300.0*rdxCp2[0]*phiUx[1]+218700.0*rdxCp2[0]*phiLy[1]+((-14029.6115413079*phiUx[0])-98207.28078915528*phiLy[0]+224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((814886.6036909672*phiUx[0]-2061183.762277152*phiLy[0]+2492594.31717237*bcVals[0])*rdxCp2Sq[0]-664200.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-2961900.0*rdxCp2R3[0]*phiUx[1])-3134700.0*rdxCp2R3[0]*phiLy[1]+(3407636.758811008*phiUx[0]-6703036.625291553*phiLy[0]+6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-31176.91453623978*rdxCp2[0]*phiUx[1])-155884.5726811989*rdxCp2[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-175759.8556980518*rdxCp2Sq[0]*phiUx[1])-687061.2540923841*rdxCp2Sq[0]*phiLy[1]+(166050.0*phiUx[0]-332100.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-244738.7791094823*rdxCp2R3[0]*phiUx[1])-469212.5637704087*rdxCp2R3[0]*phiLy[1]+(266400.0*phiUx[0]-387000.0*phiLy[0]+241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*volFac+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+1052640.0*phiLy[0]*rdxCp2R4[1]+((-893738.2167055405*rdxCp2[0]*phiUx[1])-454351.5678414679*rdxCp2[0]*phiLy[1]+(928800.0*phiUx[0]+2563956.0*phiLy[0]-1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-1772702.040022519*rdxCp2Sq[0]*phiUx[1])-108651.5471587956*rdxCp2Sq[0]*phiLy[1]+(1842246.0*phiUx[0]+532953.0*phiLy[0]-5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-439178.8027671645*rdxCp2R3[0]*phiUx[1])+1870.614872174387*rdxCp2R3[0]*phiLy[1]+(456408.0*phiUx[0]-11340.0*phiLy[0]-1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiUx[1]+(12960.0*phiUx[0]-37440.0*bcVals[0])*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-5263200.0*phiLy[1]*rdxCp2R4[1]+(6450000.0*rdxCp2[0]*phiUx[1]-1408860.0*rdxCp2[0]*phiLy[1]+((-6703036.625291553*phiUx[0])+3407636.758811008*phiLy[0]-1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(1983375.0*rdxCp2Sq[0]*phiUx[1]+63990.0*rdxCp2Sq[0]*phiLy[1]+((-2061183.762277152*phiUx[0])+814886.6036909672*phiLy[0]-1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(94500.0*rdxCp2R3[0]*phiUx[1]+9720.0*rdxCp2R3[0]*phiLy[1]+((-98207.28078915528*phiUx[0])-14029.6115413079*phiLy[0]-2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-74824.59488697546*bcVals[0]*rdxCp2R4[0]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+643491.516027989*phiLy[0]*rdxCp2R4[1]+((-154800.0*rdxCp2[0]*phiUx[1])-106560.0*rdxCp2[0]*phiLy[1]+(160872.8790069972*phiUx[0]+1745283.675738703*phiLy[0]-285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-290400.0*rdxCp2Sq[0]*phiUx[1])-66420.0*rdxCp2Sq[0]*phiLy[1]+(301792.5327108011*phiUx[0]+659547.6270141526*phiLy[0]-871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-63000.0*rdxCp2R3[0]*phiUx[1])-10800.0*rdxCp2R3[0]*phiLy[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+(372390.9236273086*rdxCp2[0]*phiUx[1]-2016730.678300898*rdxCp2[0]*phiLy[1]+((-387000.0*phiUx[0])+266400.0*phiLy[0]-688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((166050.0*phiLy[0]-580800.0*bcVals[0])*rdxCp2Sq[0]-543205.7742697513*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-25980.76211353316*rdxCp2R3[0]*phiUx[1])-18706.14872174387*rdxCp2R3[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*volFac+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-1440.0*phiLy[0]*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiUx[1]+1818.653347947321*rdxCp2[0]*phiLy[1]+((-1890.0*phiUx[0])-11799.0*phiLy[0]+3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(11353.59304361399*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLy[1]+((-11799.0*phiUx[0])-1890.0*phiLy[0]+33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiUx[1]+(4160.0*bcVals[0]-1440.0*phiUx[0])*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+1440.0*phiLy[1]*rdxCp2R3[1]+((-2625.0*rdxCp2[0]*phiUx[1])+2664.0*rdxCp2[0]*phiLy[1]+(2727.980021920981*phiUx[0]-2727.980021920981*phiLy[0]+4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiUx[1])+324.0*rdxCp2Sq[0]*phiLy[1]+(467.6537180435967*phiUx[0]-467.6537180435967*phiLy[0]+14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1662.768775266122*bcVals[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(450.0*rdxCp2[0]*phiUx[1]+450.0*rdxCp2[0]*phiLy[1]+((-467.6537180435967*phiUx[0])+467.6537180435967*phiLy[0]+831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2625.0*rdxCp2Sq[0]*phiUx[1]-450.0*rdxCp2Sq[0]*phiLy[1]+((-2727.980021920981*phiUx[0])+2727.980021920981*phiLy[0]+8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-1082.531754730548*rdxCp2[0]*phiUx[1])-4546.633369868302*rdxCp2[0]*phiLy[1]+(1125.0*phiUx[0]-1125.0*phiLy[0]+2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1082.531754730548*rdxCp2Sq[0]*phiUx[1]-779.4228634059946*rdxCp2Sq[0]*phiLy[1]+((-1125.0*phiUx[0])+1125.0*phiLy[0]+4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-328950.0*phiUy[0]*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiLx[1]-3096000.0*rdxCp2[0]*bcVals[1]+(193500.0*phiLx[0]-2715915.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiLx[1]-9081960.0*rdxCp2Sq[0]*bcVals[1]+(193500.0*phiUy[0]-2715915.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiLx[1]-1006200.0*rdxCp2R3[0]*bcVals[1]-328950.0*phiLx[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+328950.0*phiUy[1]*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiLx[1]+4468691.083527703*rdxCp2[0]*bcVals[1]+((-567939.4598018348*phiUy[0])-1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiLx[1]+4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(33515.18312645776*phiUy[0]-1919718.512568964*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiLx[1]+402182.1975174932*rdxCp2R3[0]*bcVals[1]-201091.0987587466*phiLx[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiLx[1]-536242.9300233242*rdxCp2[0]*bcVals[1]+(33515.18312645776*phiLx[0]-1919718.512568964*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiLx[1]-1098466.622160182*rdxCp2Sq[0]*bcVals[1]+((-1117172.770881926*phiUy[0])-567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiLx[1]+258000.0*rdxCp2[0]*bcVals[1]+((-44400.0*phiUy[0])-64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiLx[1]-40200.0*rdxCp2Sq[0]*bcVals[1]+((-64500.0*phiUy[0])-44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*volFac+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+12960.0*phiUy[0]*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiLx[1]+181440.0*rdxCp2[0]*bcVals[1]+(456408.0*phiUy[0]-11340.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiLx[1]+2626452.0*rdxCp2Sq[0]*bcVals[1]+(1842246.0*phiUy[0]+532953.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiLx[1]+8373744.0*rdxCp2R3[0]*bcVals[1]+(928800.0*phiUy[0]+2563956.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiLx[1]+3219840.0*rdxCp2R4[0]*bcVals[1]+1052640.0*phiLx[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*volFac+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-12960.0*phiUy[1]*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiLx[1]-261886.0821044141*rdxCp2[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiLx[1]-1922680.319449907*rdxCp2Sq[0]*bcVals[1]+(301792.5327108011*phiUy[0]+659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiLx[1]-3812313.1094914*rdxCp2R3[0]*bcVals[1]+(160872.8790069972*phiUy[0]+1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiLx[1]-1286983.032055978*rdxCp2R4[0]*bcVals[1]+643491.516027989*phiLx[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiLx[1]+224473.7846609264*rdxCp2[0]*bcVals[1]+((-98207.28078915528*phiUy[0])-14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiLx[1]+2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiLx[0]-2061183.762277152*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiLx[1]+6590799.732961089*rdxCp2R3[0]*bcVals[1]+(3407636.758811008*phiLx[0]-6703036.625291553*phiUy[0])*rdxCp2R3[0])*rdxCp2[1]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiLx[1]-108000.0*rdxCp2[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiLx[1]-332100.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiLx[1]+241200.0*rdxCp2R3[0]*bcVals[1]+(266400.0*phiLx[0]-387000.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*volFac+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1052640.0*phiUy[0]*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiLx[1]+1651200.0*rdxCp2[0]*bcVals[1]+(2563956.0*phiUy[0]+928800.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiLx[1]+5004704.0*rdxCp2Sq[0]*bcVals[1]+(532953.0*phiUy[0]+1842246.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiLx[1]+1303392.0*rdxCp2R3[0]*bcVals[1]+(456408.0*phiLx[0]-11340.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiLx[1]+37440.0*rdxCp2R4[0]*bcVals[1]+12960.0*phiLx[0]*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+5263200.0*phiUy[1]*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiLx[1]+1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(3407636.758811008*phiUy[0]-6703036.625291553*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiLx[1]+1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiUy[0]-2061183.762277152*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiLx[1]+2731097.713374604*rdxCp2R3[0]*bcVals[1]+((-14029.6115413079*phiUy[0])-98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+74824.59488697546*rdxCp2R4[0]*bcVals[1])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiLx[1]+285996.229345773*rdxCp2[0]*bcVals[1]+(1745283.675738703*phiUy[0]+160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiLx[1]+871845.0944978701*rdxCp2Sq[0]*bcVals[1]+(659547.6270141526*phiUy[0]+301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiLx[1]+201610.7140010173*rdxCp2R3[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiLx[1]+688000.0*rdxCp2[0]*bcVals[1]+(266400.0*phiUy[0]-387000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+580800.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiLx[1]+99600.0*rdxCp2R3[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*volFac+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-1440.0*phiUy[0]*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiLx[1]-3360.0*rdxCp2[0]*bcVals[1]+((-11799.0*phiUy[0])-1890.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiLx[1]-33726.0*rdxCp2Sq[0]*bcVals[1]+((-1890.0*phiUy[0])-11799.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiLx[1]-4160.0*rdxCp2R3[0]*bcVals[1]-1440.0*phiLx[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-1440.0*phiUy[1]*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiLx[1]-4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiLx[1]-14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]-1662.768775266122*rdxCp2R3[0]*bcVals[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiLx[1]-831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiLx[1]-8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiLx[1]-2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiLx[1]-4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiUy[0]-1125.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+328950.0*phiLy[0]*rdxCp2R3[1]+((-1533436.541464953*rdxCp2[0]*phiLy[1])-335151.8312645776*rdxCp2[0]*phiLx[1]+3096000.0*rdxCp2[0]*bcVals[1]+(2715915.0*phiLy[0]-193500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(90490.99444143593*rdxCp2Sq[0]*phiLy[1]+2176434.423012786*rdxCp2Sq[0]*phiLx[1]+9081960.0*rdxCp2Sq[0]*bcVals[1]+(2715915.0*phiLx[0]-193500.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiLx[1]+1006200.0*rdxCp2R3[0]*bcVals[1]+328950.0*phiLx[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]-328950.0*phiLy[1]*rdxCp2R3[1]+((-125505.0*rdxCp2[0]*phiLy[1])+1290000.0*rdxCp2[0]*phiLx[1]-4468691.083527703*rdxCp2[0]*bcVals[1]+(567939.4598018348*phiLy[0]+1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(40635.0*rdxCp2Sq[0]*phiLy[1]+2054550.0*rdxCp2Sq[0]*phiLx[1]-4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(1919718.512568964*phiLx[0]-33515.18312645776*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiLx[1]-402182.1975174932*rdxCp2R3[0]*bcVals[1]+201091.0987587466*phiLx[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-359640.0*rdxCp2[0]*phiLy[1])-58050.0*rdxCp2[0]*phiLx[1]+536242.9300233242*rdxCp2[0]*bcVals[1]+(1919718.512568964*phiLy[0]-33515.18312645776*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-522450.0*rdxCp2Sq[0]*phiLy[1])+493650.0*rdxCp2Sq[0]*phiLx[1]+1098466.622160182*rdxCp2Sq[0]*bcVals[1]+(1117172.770881926*phiLy[0]+567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+((-437394.7904353685*rdxCp2[0]*phiLy[1])+74478.18472546172*rdxCp2[0]*phiLx[1]-258000.0*rdxCp2[0]*bcVals[1]+(44400.0*phiLy[0]+64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-78202.0939617348*rdxCp2Sq[0]*phiLy[1])+40789.79651824706*rdxCp2Sq[0]*phiLx[1]+40200.0*rdxCp2Sq[0]*bcVals[1]+(64500.0*phiLy[0]+44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*volFac+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-12960.0*phiLy[0]*rdxCp2R4[1]+(176773.1054204795*rdxCp2[0]*phiLy[1]+19641.45615783106*rdxCp2[0]*phiLx[1]-181440.0*rdxCp2[0]*bcVals[1]+(11340.0*phiLx[0]-456408.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+(814839.8383191626*rdxCp2Sq[0]*phiLy[1]-386469.0325912283*rdxCp2Sq[0]*phiLx[1]-2626452.0*rdxCp2Sq[0]*bcVals[1]+((-1842246.0*phiLy[0])-532953.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(434356.7733188925*rdxCp2R3[0]*phiLy[1]-2064285.865273508*rdxCp2R3[0]*phiLx[1]-8373744.0*rdxCp2R3[0]*bcVals[1]+((-928800.0*phiLy[0])-2563956.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiLx[1]-3219840.0*rdxCp2R4[0]*bcVals[1]-1052640.0*phiLx[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*volFac+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+12960.0*phiLy[1]*rdxCp2R4[1]+(157788.0*rdxCp2[0]*phiLy[1]-75600.0*rdxCp2[0]*phiLx[1]+261886.0821044141*rdxCp2[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(465750.0*rdxCp2Sq[0]*phiLy[1]-727155.0*rdxCp2Sq[0]*phiLx[1]+1922680.319449907*rdxCp2Sq[0]*bcVals[1]+((-301792.5327108011*phiLy[0])-659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(195048.0*rdxCp2R3[0]*phiLy[1]-1862820.0*rdxCp2R3[0]*phiLx[1]+3812313.1094914*rdxCp2R3[0]*bcVals[1]+((-160872.8790069972*phiLy[0])-1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiLx[1]+1286983.032055978*rdxCp2R4[0]*bcVals[1]-643491.516027989*phiLx[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+(218700.0*rdxCp2[0]*phiLy[1]+24300.0*rdxCp2[0]*phiLx[1]-224473.7846609264*rdxCp2[0]*bcVals[1]+(98207.28078915528*phiLy[0]+14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-664200.0*rdxCp2Sq[0]*phiLx[1])-2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLy[0]-814886.6036909672*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3134700.0*rdxCp2R3[0]*phiLy[1])-2961900.0*rdxCp2R3[0]*phiLx[1]-6590799.732961089*rdxCp2R3[0]*bcVals[1]+(6703036.625291553*phiLy[0]-3407636.758811008*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-155884.5726811989*rdxCp2[0]*phiLy[1])-31176.91453623978*rdxCp2[0]*phiLx[1]+108000.0*rdxCp2[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-687061.2540923841*rdxCp2Sq[0]*phiLy[1])-175759.8556980518*rdxCp2Sq[0]*phiLx[1]+332100.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-469212.5637704087*rdxCp2R3[0]*phiLy[1])-244738.7791094823*rdxCp2R3[0]*phiLx[1]-241200.0*rdxCp2R3[0]*bcVals[1]+(387000.0*phiLy[0]-266400.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*volFac+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1052640.0*phiLy[0]*rdxCp2R4[1]+((-454351.5678414679*rdxCp2[0]*phiLy[1])-893738.2167055405*rdxCp2[0]*phiLx[1]-1651200.0*rdxCp2[0]*bcVals[1]+((-2563956.0*phiLy[0])-928800.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-108651.5471587956*rdxCp2Sq[0]*phiLy[1])-1772702.040022519*rdxCp2Sq[0]*phiLx[1]-5004704.0*rdxCp2Sq[0]*bcVals[1]+((-532953.0*phiLy[0])-1842246.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(1870.614872174387*rdxCp2R3[0]*phiLy[1]-439178.8027671645*rdxCp2R3[0]*phiLx[1]-1303392.0*rdxCp2R3[0]*bcVals[1]+(11340.0*phiLy[0]-456408.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiLx[1]-37440.0*rdxCp2R4[0]*bcVals[1]-12960.0*phiLx[0]*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]-5263200.0*phiLy[1]*rdxCp2R4[1]+((-1408860.0*rdxCp2[0]*phiLy[1])+6450000.0*rdxCp2[0]*phiLx[1]-1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(6703036.625291553*phiLx[0]-3407636.758811008*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+(63990.0*rdxCp2Sq[0]*phiLy[1]+1983375.0*rdxCp2Sq[0]*phiLx[1]-1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLx[0]-814886.6036909672*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(9720.0*rdxCp2R3[0]*phiLy[1]+94500.0*rdxCp2R3[0]*phiLx[1]-2731097.713374604*rdxCp2R3[0]*bcVals[1]+(14029.6115413079*phiLy[0]+98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-74824.59488697546*rdxCp2R4[0]*bcVals[1]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]-643491.516027989*phiLy[0]*rdxCp2R4[1]+((-106560.0*rdxCp2[0]*phiLy[1])-154800.0*rdxCp2[0]*phiLx[1]-285996.229345773*rdxCp2[0]*bcVals[1]+((-1745283.675738703*phiLy[0])-160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-66420.0*rdxCp2Sq[0]*phiLy[1])-290400.0*rdxCp2Sq[0]*phiLx[1]-871845.0944978701*rdxCp2Sq[0]*bcVals[1]+((-659547.6270141526*phiLy[0])-301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-10800.0*rdxCp2R3[0]*phiLy[1])-63000.0*rdxCp2R3[0]*phiLx[1]-201610.7140010173*rdxCp2R3[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+((-2016730.678300898*rdxCp2[0]*phiLy[1])+372390.9236273086*rdxCp2[0]*phiLx[1]-688000.0*rdxCp2[0]*bcVals[1]+(387000.0*phiLx[0]-266400.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-543205.7742697513*rdxCp2Sq[0]*phiLy[1])-580800.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-18706.14872174387*rdxCp2R3[0]*phiLy[1])-25980.76211353316*rdxCp2R3[0]*phiLx[1]-99600.0*rdxCp2R3[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*volFac+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+1440.0*phiLy[0]*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiLy[1]+1818.653347947321*rdxCp2[0]*phiLx[1]+3360.0*rdxCp2[0]*bcVals[1]+(11799.0*phiLy[0]+1890.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2Sq[0]*phiLy[1]+11353.59304361399*rdxCp2Sq[0]*phiLx[1]+33726.0*rdxCp2Sq[0]*bcVals[1]+(1890.0*phiLy[0]+11799.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiLx[1]+4160.0*rdxCp2R3[0]*bcVals[1]+1440.0*phiLx[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+1440.0*phiLy[1]*rdxCp2R3[1]+(2664.0*rdxCp2[0]*phiLy[1]-2625.0*rdxCp2[0]*phiLx[1]+4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLy[0]-2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(324.0*rdxCp2Sq[0]*phiLy[1]-450.0*rdxCp2Sq[0]*phiLx[1]+14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLy[0]-467.6537180435967*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+1662.768775266122*rdxCp2R3[0]*bcVals[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+(450.0*rdxCp2[0]*phiLy[1]+450.0*rdxCp2[0]*phiLx[1]+831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiLy[1])+2625.0*rdxCp2Sq[0]*phiLx[1]+8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-4546.633369868302*rdxCp2[0]*phiLy[1])-1082.531754730548*rdxCp2[0]*phiLx[1]+2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLy[0]-1125.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-779.4228634059946*rdxCp2Sq[0]*phiLy[1])+1082.531754730548*rdxCp2Sq[0]*phiLx[1]+4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 
  double *dxUy = dx[3]; 
  double *dxLy = dx[4]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdxCp2[2]; 
  double rdxLx[2]; 
  double rdxUx[2]; 
  double rdxLxSq[2]; 
  double rdxUxSq[2]; 
  double rdxLxCu[2]; 
  double rdxUxCu[2]; 
  double rdxLxR4[2]; 
  double rdxUxR4[2]; 
  double rdxLy[2]; 
  double rdxUy[2]; 
  double rdxLySq[2]; 
  double rdxUySq[2]; 
  double rdxLyCu[2]; 
  double rdxUyCu[2]; 
  double rdxLyR4[2]; 
  double rdxUyR4[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = volFac*4.0/(dxLx[1]*dxLx[1]); 
  rdxUx[1]   = volFac*4.0/(dxUx[1]*dxUx[1]); 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = volFac*4.0/(dxLy[0]*dxLy[0]); 
  rdxUy[0]   = volFac*4.0/(dxUy[0]*dxUy[0]); 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = volFac*4.0/(dxLy[1]*dxLy[1]); 
  rdxUy[1]   = volFac*4.0/(dxUy[1]*dxUy[1]); 
  rdxLySq[1] = rdxLy[1]*rdxLy[1]; 
  rdxUySq[1] = rdxUy[1]*rdxUy[1]; 
  rdxLyCu[1] = rdxLySq[1]*rdxLy[1]; 
  rdxUyCu[1] = rdxUySq[1]*rdxUy[1]; 
  rdxLyR4[1] = rdxLyCu[1]*rdxLy[1]; 
  rdxUyR4[1] = rdxUyCu[1]*rdxUy[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((((9600.0*rdxUx[0]-9600.0*rdxLx[0])*rdxUySq[1]+(9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUy[1]+(9600.0*rdxLx[0]-9600.0*rdxUx[0])*rdxLySq[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLy[1])*rho[3]+((-415.6921938165305*rdxUyCu[1])+((-27435.68479189101*rdxLy[1])-25495.78788741387*rdxUx[0]-25495.78788741387*rdxLx[0])*rdxUySq[1]+(27435.68479189101*rdxLySq[1]-25080.09569359734*rdxUxSq[0]-23140.1987891202*rdxLx[0]*rdxUx[0]-25080.09569359734*rdxLxSq[0])*rdxUy[1]+415.6921938165305*rdxLyCu[1]+(25495.78788741387*rdxUx[0]+25495.78788741387*rdxLx[0])*rdxLySq[1]+(25080.09569359734*rdxUxSq[0]+23140.1987891202*rdxLx[0]*rdxUx[0]+25080.09569359734*rdxLxSq[0])*rdxLy[1])*rho[2]+((25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxUySq[1]+((23140.1987891202*rdxLx[0]-23140.1987891202*rdxUx[0])*rdxLy[1]-25495.78788741387*rdxUxSq[0]+25495.78788741387*rdxLxSq[0])*rdxUy[1]+(25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxLySq[1]+(25495.78788741387*rdxLxSq[0]-25495.78788741387*rdxUxSq[0])*rdxLy[1]-415.6921938165305*rdxUxCu[0]-27435.68479189101*rdxLx[0]*rdxUxSq[0]+27435.68479189101*rdxLxSq[0]*rdxUx[0]+415.6921938165305*rdxLxCu[0])*rho[1]+1104.0*rho[0]*rdxUyCu[1]+(75072.0*rho[0]*rdxLy[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0])*rdxUySq[1]+(75072.0*rho[0]*rdxLySq[1]+(164368.0*rdxUx[0]+164368.0*rdxLx[0])*rho[0]*rdxLy[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0])*rdxUy[1]+1104.0*rho[0]*rdxLyCu[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0]*rdxLySq[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0]*rdxLy[1]+(1104.0*rdxUxCu[0]+75072.0*rdxLx[0]*rdxUxSq[0]+75072.0*rdxLxSq[0]*rdxUx[0]+1104.0*rdxLxCu[0])*rho[0])*omega*volFac+(((9375.0*rdxUx[0]-9375.0*rdxLx[0])*rdxUyCu[1]+((12525.0*rdxUx[0]-12525.0*rdxLx[0])*rdxLy[1]+9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+((17775.0*rdxUx[0]-17775.0*rdxLx[0])*rdxLySq[1]+(18000.0*rdxUxSq[0]-18000.0*rdxLxSq[0])*rdxLy[1]+225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]+9600.0*rdxUxSq[0]+18000.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])+9375.0*rdxUxCu[0]+12525.0*rdxLx[0]*rdxUxSq[0]+17775.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+((-9600.0*rdxUxSq[0])-18000.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-9375.0*rdxUxCu[0])-12525.0*rdxLx[0]*rdxUxSq[0]-17775.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((17775.0*rdxLx[0]-17775.0*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((12525.0*rdxLx[0]-12525.0*rdxUx[0])*rdxLySq[1]+(18000.0*rdxLxSq[0]-18000.0*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9375.0*rdxLx[0]-9375.0*rdxUx[0])*rdxLyCu[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLySq[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-18000.0*rdxLx[0]*rdxUx[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-17775.0*rdxLx[0]*rdxUxSq[0]-12525.0*rdxLxSq[0]*rdxUx[0]-9375.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(18000.0*rdxLx[0]*rdxUx[0]+9600.0*rdxLxSq[0])*rdxLySq[1]+(17775.0*rdxLx[0]*rdxUxSq[0]+12525.0*rdxLxSq[0]*rdxUx[0]+9375.0*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((-415.6921938165305*rdxUyR4[1])+((-28630.79984911354*rdxLy[1])-25729.61474643567*rdxUx[0]-25729.61474643567*rdxLx[0])*rdxUyCu[1]+((-52637.02404201817*rdxLySq[1])+((-88966.78973077537*rdxUx[0])-88966.78973077537*rdxLx[0])*rdxLy[1]-25911.4800812304*rdxUxSq[0]-78842.95276053529*rdxLx[0]*rdxUx[0]-25911.4800812304*rdxLxSq[0])*rdxUySq[1]+((-779.4228634059946*rdxLyCu[1])+((-48038.4291479228*rdxUx[0])-48038.4291479228*rdxLx[0])*rdxLySq[1]+((-47856.56381312806*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0]-47856.56381312806*rdxLxSq[0])*rdxLy[1]-597.5575286112626*rdxUxCu[0]-40633.91194556586*rdxLx[0]*rdxUxSq[0]-40633.91194556586*rdxLxSq[0]*rdxUx[0]-597.5575286112626*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-233.8268590217983*rdxUx[0]*rdxUyCu[1])+((-15432.57269543869*rdxUx[0]*rdxLy[1])-9145.22826396367*rdxUxSq[0]-19537.53310937693*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15432.57269543869*rdxUx[0]*rdxLySq[1]-8911.401404941873*rdxUxCu[0]-13016.36181888011*rdxLx[0]*rdxUxSq[0]-19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]+19537.53310937693*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8911.401404941873*rdxUxCu[0]+13016.36181888011*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+(779.4228634059946*rdxLy[1]*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(48038.4291479228*rdxUx[0]+48038.4291479228*rdxLx[0])*rdxLy[1])*rdxUySq[1]+(28630.79984911354*rdxLyCu[1]+(88966.78973077537*rdxUx[0]+88966.78973077537*rdxLx[0])*rdxLySq[1]+(47856.56381312806*rdxUxSq[0]+99090.62670101546*rdxLx[0]*rdxUx[0]+47856.56381312806*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+415.6921938165305*rdxLyR4[1]+(25729.61474643567*rdxUx[0]+25729.61474643567*rdxLx[0])*rdxLyCu[1]+(25911.4800812304*rdxUxSq[0]+78842.95276053529*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*rdxLySq[1]+(597.5575286112626*rdxUxCu[0]+40633.91194556586*rdxLx[0]*rdxUxSq[0]+40633.91194556586*rdxLxSq[0]*rdxUx[0]+597.5575286112626*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15432.57269543869*rdxLx[0]*rdxLy[1])-19537.53310937693*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(15432.57269543869*rdxLx[0]*rdxLySq[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-13016.36181888011*rdxLxSq[0]*rdxUx[0]-8911.401404941873*rdxLxCu[0])*rdxUy[1]+233.8268590217983*rdxLx[0]*rdxLyCu[1]+(19537.53310937693*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(19303.70625035514*rdxLx[0]*rdxUxSq[0]+13016.36181888011*rdxLxSq[0]*rdxUx[0]+8911.401404941873*rdxLxCu[0])*rdxLy[1])*phiLx[2]+(396.0*phiUy[0]-36.0*phiC[0])*rdxUyR4[1]+((27378.0*phiUy[0]+846.0*phiLy[0]-4824.0*phiC[0])*rdxLy[1]+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiUy[1]-597.5575286112626*rdxUx[0]*phiUx[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(24531.0*phiUy[0]+621.0*phiUx[0]-3072.0*phiC[0])*rdxUx[0]+(24531.0*phiUy[0]+621.0*phiLx[0]-3072.0*phiC[0])*rdxLx[0])*rdxUyCu[1]+((57078.0*phiUy[0]+57078.0*phiLy[0]-161676.0*phiC[0])*rdxLySq[1]+((13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(92457.0*phiUy[0]+42228.0*phiUx[0]+52131.0*phiLy[0]-208896.0*phiC[0])*rdxUx[0]+(92457.0*phiUy[0]+52131.0*phiLy[0]+42228.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0])*rdxLy[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiUy[1]+((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUy[0]+24756.0*phiUx[0]-6072.0*phiC[0])*rdxUxSq[0]+(79932.0*phiUy[0]+51906.0*phiUx[0]+51906.0*phiLx[0]-207144.0*phiC[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiUy[0]+24756.0*phiLx[0]-6072.0*phiC[0])*rdxLxSq[0])*rdxUySq[1]+((846.0*phiUy[0]+27378.0*phiLy[0]-4824.0*phiC[0])*rdxLyCu[1]+((19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(52131.0*phiUy[0]+42228.0*phiUx[0]+92457.0*phiLy[0]-208896.0*phiC[0])*rdxUx[0]+(52131.0*phiUy[0]+92457.0*phiLy[0]+42228.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0])*rdxLySq[1]+((19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiUy[1]+((-78842.95276053529*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0])*phiUx[1]+(19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiLy[1]+(99090.62670101546*rdxLx[0]*rdxUx[0]+78842.95276053529*rdxLxSq[0])*phiLx[1]+(51906.0*phiUy[0]+79932.0*phiUx[0]+51906.0*phiLy[0]-207144.0*phiC[0])*rdxUxSq[0]+(104982.0*phiUy[0]+104982.0*phiUx[0]+104982.0*phiLy[0]+104982.0*phiLx[0]-500088.0*phiC[0])*rdxLx[0]*rdxUx[0]+(51906.0*phiUy[0]+51906.0*phiLy[0]+79932.0*phiLx[0]-207144.0*phiC[0])*rdxLxSq[0])*rdxLy[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiUy[1]+((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(621.0*phiUy[0]+24531.0*phiUx[0]-3072.0*phiC[0])*rdxUxCu[0]+(42228.0*phiUy[0]+92457.0*phiUx[0]+52131.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0]*rdxUxSq[0]+(42228.0*phiUy[0]+52131.0*phiUx[0]+92457.0*phiLx[0]-208896.0*phiC[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiUy[0]+24531.0*phiLx[0]-3072.0*phiC[0])*rdxLxCu[0])*rdxUy[1]+(396.0*phiLy[0]-36.0*phiC[0])*rdxLyR4[1]+((-597.5575286112626*rdxUx[0]*phiUx[1])+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiLy[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(621.0*phiUx[0]+24531.0*phiLy[0]-3072.0*phiC[0])*rdxUx[0]+(24531.0*phiLy[0]+621.0*phiLx[0]-3072.0*phiC[0])*rdxLx[0])*rdxLyCu[1]+(((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiLy[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUx[0]+24756.0*phiLy[0]-6072.0*phiC[0])*rdxUxSq[0]+(51906.0*phiUx[0]+79932.0*phiLy[0]+51906.0*phiLx[0]-207144.0*phiC[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiLy[0]+24756.0*phiLx[0]-6072.0*phiC[0])*rdxLxSq[0])*rdxLySq[1]+(((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiLy[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(24531.0*phiUx[0]+621.0*phiLy[0]-3072.0*phiC[0])*rdxUxCu[0]+(92457.0*phiUx[0]+42228.0*phiLy[0]+52131.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0]*rdxUxSq[0]+(52131.0*phiUx[0]+42228.0*phiLy[0]+92457.0*phiLx[0]-208896.0*phiC[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiLy[0]+24531.0*phiLx[0]-3072.0*phiC[0])*rdxLxCu[0])*rdxLy[1]+((-415.6921938165305*rdxUxR4[0])-28630.79984911354*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(779.4228634059946*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+28630.79984911354*rdxLxCu[0]*rdxUx[0]+415.6921938165305*rdxLxR4[0])*phiLx[1]+(396.0*phiUx[0]-36.0*phiC[0])*rdxUxR4[0]+(27378.0*phiUx[0]+846.0*phiLx[0]-4824.0*phiC[0])*rdxLx[0]*rdxUxCu[0]+(57078.0*phiUx[0]+57078.0*phiLx[0]-161676.0*phiC[0])*rdxLxSq[0]*rdxUxSq[0]+(846.0*phiUx[0]+27378.0*phiLx[0]-4824.0*phiC[0])*rdxLxCu[0]*rdxUx[0]+(396.0*phiLx[0]-36.0*phiC[0])*rdxLxR4[0])*omega+36.0*phiC[0]*rdxUyR4[1]+(4824.0*phiC[0]*rdxLy[1]+3072.0*phiC[0]*rdxUx[0]+3072.0*phiC[0]*rdxLx[0])*rdxUyCu[1]+(161676.0*phiC[0]*rdxLySq[1]+(208896.0*phiC[0]*rdxUx[0]+208896.0*phiC[0]*rdxLx[0])*rdxLy[1]+6072.0*phiC[0]*rdxUxSq[0]+207144.0*phiC[0]*rdxLx[0]*rdxUx[0]+6072.0*phiC[0]*rdxLxSq[0])*rdxUySq[1]+(4824.0*phiC[0]*rdxLyCu[1]+(208896.0*phiC[0]*rdxUx[0]+208896.0*phiC[0]*rdxLx[0])*rdxLySq[1]+(207144.0*phiC[0]*rdxUxSq[0]+500088.0*phiC[0]*rdxLx[0]*rdxUx[0]+207144.0*phiC[0]*rdxLxSq[0])*rdxLy[1]+3072.0*phiC[0]*rdxUxCu[0]+208896.0*phiC[0]*rdxLx[0]*rdxUxSq[0]+208896.0*phiC[0]*rdxLxSq[0]*rdxUx[0]+3072.0*phiC[0]*rdxLxCu[0])*rdxUy[1]+36.0*phiC[0]*rdxLyR4[1]+(3072.0*phiC[0]*rdxUx[0]+3072.0*phiC[0]*rdxLx[0])*rdxLyCu[1]+(6072.0*phiC[0]*rdxUxSq[0]+207144.0*phiC[0]*rdxLx[0]*rdxUx[0]+6072.0*phiC[0]*rdxLxSq[0])*rdxLySq[1]+(3072.0*phiC[0]*rdxUxCu[0]+208896.0*phiC[0]*rdxLx[0]*rdxUxSq[0]+208896.0*phiC[0]*rdxLxSq[0]*rdxUx[0]+3072.0*phiC[0]*rdxLxCu[0])*rdxLy[1]+36.0*phiC[0]*rdxUxR4[0]+4824.0*phiC[0]*rdxLx[0]*rdxUxCu[0]+161676.0*phiC[0]*rdxLxSq[0]*rdxUxSq[0]+4824.0*phiC[0]*rdxLxCu[0]*rdxUx[0]+36.0*phiC[0]*rdxLxR4[0])/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxUyCu[1]+(27435.68479189101*rdxLy[1]+9976.61265159673*rdxUx[0]+9976.61265159673*rdxLx[0])*rdxUySq[1]+((-27435.68479189101*rdxLySq[1])+9560.920457780201*rdxUxSq[0]-7898.15168251408*rdxLx[0]*rdxUx[0]+9560.920457780201*rdxLxSq[0])*rdxUy[1]-415.6921938165305*rdxLyCu[1]+((-9976.61265159673*rdxUx[0])-9976.61265159673*rdxLx[0])*rdxLySq[1]+((-9560.920457780201*rdxUxSq[0])+7898.15168251408*rdxLx[0]*rdxUx[0]-9560.920457780201*rdxLxSq[0])*rdxLy[1])*rho[3]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[2]+((-1104.0*rdxUyCu[1])+((-75072.0*rdxLy[1])-27600.0*rdxUx[0]-27600.0*rdxLx[0])*rdxUySq[1]+((-75072.0*rdxLySq[1])+((-126960.0*rdxUx[0])-126960.0*rdxLx[0])*rdxLy[1]-26928.0*rdxUxSq[0]-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxUy[1]-1104.0*rdxLyCu[1]+((-27600.0*rdxUx[0])-27600.0*rdxLx[0])*rdxLySq[1]+((-26928.0*rdxUxSq[0])-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxLy[1]-432.0*rdxUxCu[0]-29376.0*rdxLx[0]*rdxUxSq[0]-29376.0*rdxLxSq[0]*rdxUx[0]-432.0*rdxLxCu[0])*rho[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxUySq[1]+((60164.51685171252*rdxUx[0]-60164.51685171252*rdxLx[0])*rho[0]*rdxLy[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0])*rdxUy[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxLySq[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0]*rdxLy[1]+(1080.799703922979*rdxUxCu[0]+71332.78045891662*rdxLx[0]*rdxUxSq[0]-71332.78045891662*rdxLxSq[0]*rdxUx[0]-1080.799703922979*rdxLxCu[0])*rho[0])*omega*volFac+((415.6921938165305*rdxUyR4[1]+(28630.79984911354*rdxLy[1]+10574.170180208*rdxUx[0]+10574.170180208*rdxLx[0])*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(68719.1157902952*rdxUx[0]+68719.1157902952*rdxLx[0])*rdxLy[1]+10392.30484541326*rdxUxSq[0]+47804.60228890101*rdxLx[0]*rdxUx[0]+10392.30484541326*rdxLxSq[0])*rdxUySq[1]+(779.4228634059946*rdxLyCu[1]+(19303.70625035514*rdxUx[0]+19303.70625035514*rdxLx[0])*rdxLySq[1]+(18758.11024597094*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0]+18758.11024597094*rdxLxSq[0])*rdxLy[1]+233.8268590217983*rdxUxCu[0]+15900.22641348229*rdxLx[0]*rdxUxSq[0]+15900.22641348229*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-181.8653347947321*rdxUx[0]*rdxUyCu[1])+((-12003.11209645232*rdxUx[0]*rdxLy[1])+9145.22826396367*rdxUxSq[0]-17874.76433411081*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(12003.11209645232*rdxUx[0]*rdxLySq[1]+9327.093598758403*rdxUxCu[0]+3455.441361099909*rdxLx[0]*rdxUxSq[0]-17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+181.8653347947321*rdxUx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxUxSq[0])*rdxLySq[1]+((-9327.093598758403*rdxUxCu[0])-3455.441361099909*rdxLx[0]*rdxUxSq[0]+17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((-779.4228634059946*rdxLy[1]*rdxUyCu[1])+(((-19303.70625035514*rdxUx[0])-19303.70625035514*rdxLx[0])*rdxLy[1]-52637.02404201817*rdxLySq[1])*rdxUySq[1]+((-28630.79984911354*rdxLyCu[1])+((-68719.1157902952*rdxUx[0])-68719.1157902952*rdxLx[0])*rdxLySq[1]+((-18758.11024597094*rdxUxSq[0])-40893.71956670118*rdxLx[0]*rdxUx[0]-18758.11024597094*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-415.6921938165305*rdxLyR4[1]+((-10574.170180208*rdxUx[0])-10574.170180208*rdxLx[0])*rdxLyCu[1]+((-10392.30484541326*rdxUxSq[0])-47804.60228890101*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-233.8268590217983*rdxUxCu[0])-15900.22641348229*rdxLx[0]*rdxUxSq[0]-15900.22641348229*rdxLxSq[0]*rdxUx[0]-233.8268590217983*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-181.8653347947321*rdxLx[0]*rdxUyCu[1])+((-12003.11209645232*rdxLx[0]*rdxLy[1])-17874.76433411081*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(12003.11209645232*rdxLx[0]*rdxLySq[1]-17692.89899931608*rdxLx[0]*rdxUxSq[0]+3455.441361099909*rdxLxSq[0]*rdxUx[0]+9327.093598758403*rdxLxCu[0])*rdxUy[1]+181.8653347947321*rdxLx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(17692.89899931608*rdxLx[0]*rdxUxSq[0]-3455.441361099909*rdxLxSq[0]*rdxUx[0]-9327.093598758403*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((24375.0*rdxLx[0]-24375.0*rdxUx[0])*rdxUyCu[1]+((32565.0*rdxLx[0]-32565.0*rdxUx[0])*rdxLy[1]-24960.0*rdxUxSq[0]+24960.0*rdxLxSq[0])*rdxUySq[1]+((46215.0*rdxLx[0]-46215.0*rdxUx[0])*rdxLySq[1]+(46800.0*rdxLxSq[0]-46800.0*rdxUxSq[0])*rdxLy[1]-585.0*rdxUxCu[0]-38610.0*rdxLx[0]*rdxUxSq[0]+38610.0*rdxLxSq[0]*rdxUx[0]+585.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]-8640.0*rdxUxSq[0]+19440.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])-8865.0*rdxUxCu[0]-4275.0*rdxLx[0]*rdxUxSq[0]+19215.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+(8640.0*rdxUxSq[0]-19440.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8865.0*rdxUxCu[0]+4275.0*rdxLx[0]*rdxUxSq[0]-19215.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+((46215.0*rdxUx[0]-46215.0*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((32565.0*rdxUx[0]-32565.0*rdxLx[0])*rdxLySq[1]+(46800.0*rdxUxSq[0]-46800.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(24375.0*rdxUx[0]-24375.0*rdxLx[0])*rdxLyCu[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLySq[1]+(585.0*rdxUxCu[0]+38610.0*rdxLx[0]*rdxUxSq[0]-38610.0*rdxLxSq[0]*rdxUx[0]-585.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-19440.0*rdxLx[0]*rdxUx[0]+8640.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-19215.0*rdxLx[0]*rdxUxSq[0]+4275.0*rdxLxSq[0]*rdxUx[0]+8865.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(19440.0*rdxLx[0]*rdxUx[0]-8640.0*rdxLxSq[0])*rdxLySq[1]+(19215.0*rdxLx[0]*rdxUxSq[0]-4275.0*rdxLxSq[0]*rdxUx[0]-8865.0*rdxLxCu[0])*rdxLy[1])*phiLx[2]+(36.0*phiC[1]-396.0*phiUy[1])*rdxUyR4[1]+(((-27378.0*phiUy[1])-846.0*phiLy[1]+4824.0*phiC[1])*rdxLy[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiUy[1]+483.0*rdxUx[0]*phiUx[1]+483.0*rdxLx[0]*phiLx[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*phiC[1]+(23169.64365284887*phiUy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiUy[0])*rdxLx[0])*rdxUyCu[1]+(((-57078.0*phiUy[1])-57078.0*phiLy[1]+161676.0*phiC[1])*rdxLySq[1]+(((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*phiC[1]+(33842.54072908828*phiUy[0]-40633.91194556586*phiUx[0]+50189.63625092335*phiLy[0])*rdxUx[0]+((-33842.54072908828*phiUy[0])-50189.63625092335*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLy[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiUy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*phiC[1]+(23777.59348630554*phiUy[0]+21740.70173660454*phiUx[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiUy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-846.0*phiUy[1])-27378.0*phiLy[1]+4824.0*phiC[1])*rdxLyCu[1]+(((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*phiC[1]+(50189.63625092335*phiUy[0]-40633.91194556586*phiUx[0]+33842.54072908828*phiLy[0])*rdxUx[0]+((-50189.63625092335*phiUy[0])-33842.54072908828*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLySq[1]+(((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiUy[1]+(22980.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0])*phiUx[1]+((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiLy[1]+(88110.0*rdxLx[0]*rdxUx[0]+22980.0*rdxLxSq[0])*phiLx[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*phiC[1]+(50797.58608438003*phiUy[0]-34876.57506120691*phiUx[0]+50797.58608438003*phiLy[0])*rdxUxSq[0]+(102561.6565193835*phiLx[0]-102561.6565193835*phiUx[0])*rdxLx[0]*rdxUx[0]+((-50797.58608438003*phiUy[0])-50797.58608438003*phiLy[0]+34876.57506120691*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiUy[1]+((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*phiC[1]+(607.9498334566757*phiUy[0]+22712.38223965068*phiUx[0])*rdxUxCu[0]+(40124.68900814059*phiUy[0]-44349.16092780109*phiUx[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40124.68900814059*phiUy[0])-51862.79733103487*phiUx[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiUy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxUy[1]+(36.0*phiC[1]-396.0*phiLy[1])*rdxLyR4[1]+(483.0*rdxUx[0]*phiUx[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiLy[1]+483.0*rdxLx[0]*phiLx[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*phiC[1]+(23169.64365284887*phiLy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiLy[0])*rdxLx[0])*rdxLyCu[1]+((47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiLy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*phiC[1]+(21740.70173660454*phiUx[0]+23777.59348630554*phiLy[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiLy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiLy[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*phiC[1]+(22712.38223965068*phiUx[0]+607.9498334566757*phiLy[0])*rdxUxCu[0]+((-44349.16092780109*phiUx[0])+40124.68900814059*phiLy[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-51862.79733103487*phiUx[0])-40124.68900814059*phiLy[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiLy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-396.0*rdxUxR4[0])-25758.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]+774.0*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(774.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]-25758.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[1]+(36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0])*phiC[1]+374.1229744348773*phiUx[0]*rdxUxR4[0]+(24224.46259465831*phiUx[0]+841.7766924784738*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(56024.91542162288*phiLx[0]-56024.91542162288*phiUx[0])*rdxLxSq[0]*rdxUxSq[0]+((-841.7766924784738*phiUx[0])-24224.46259465831*phiLx[0])*rdxLxCu[0]*rdxUx[0]-374.1229744348773*phiLx[0]*rdxLxR4[0])*omega-36.0*phiC[1]*rdxUyR4[1]+(((-3072.0*rdxUx[0])-3072.0*rdxLx[0])*phiC[1]-4824.0*phiC[1]*rdxLy[1])*rdxUyCu[1]+((-161676.0*phiC[1]*rdxLySq[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*phiC[1]*rdxLy[1]+((-6072.0*rdxUxSq[0])-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*phiC[1])*rdxUySq[1]+((-4824.0*phiC[1]*rdxLyCu[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*phiC[1]*rdxLySq[1]+((-207144.0*rdxUxSq[0])-500088.0*rdxLx[0]*rdxUx[0]-207144.0*rdxLxSq[0])*phiC[1]*rdxLy[1]+((-3072.0*rdxUxCu[0])-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*phiC[1])*rdxUy[1]-36.0*phiC[1]*rdxLyR4[1]+((-3072.0*rdxUx[0])-3072.0*rdxLx[0])*phiC[1]*rdxLyCu[1]+((-6072.0*rdxUxSq[0])-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*phiC[1]*rdxLySq[1]+((-3072.0*rdxUxCu[0])-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*phiC[1]*rdxLy[1]+((-36.0*rdxUxR4[0])-4824.0*rdxLx[0]*rdxUxCu[0]-161676.0*rdxLxSq[0]*rdxUxSq[0]-4824.0*rdxLxCu[0]*rdxUx[0]-36.0*rdxLxR4[0])*phiC[1]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[2] = -(1.0*((((9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxUySq[1]+((7898.15168251408*rdxLx[0]-7898.15168251408*rdxUx[0])*rdxLy[1]+9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxUy[1]+(9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxLySq[1]+(9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxLy[1]+415.6921938165305*rdxUxCu[0]+27435.68479189101*rdxLx[0]*rdxUxSq[0]-27435.68479189101*rdxLxSq[0]*rdxUx[0]-415.6921938165305*rdxLxCu[0])*rho[3]+((-432.0*rdxUyCu[1])+((-29376.0*rdxLy[1])-26928.0*rdxUx[0]-26928.0*rdxLx[0])*rdxUySq[1]+((-29376.0*rdxLySq[1])+((-81936.0*rdxUx[0])-81936.0*rdxLx[0])*rdxLy[1]-27600.0*rdxUxSq[0]-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxUy[1]-432.0*rdxLyCu[1]+((-26928.0*rdxUx[0])-26928.0*rdxLx[0])*rdxLySq[1]+((-27600.0*rdxUxSq[0])-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxLy[1]-1104.0*rdxUxCu[0]-75072.0*rdxLx[0]*rdxUxSq[0]-75072.0*rdxLxSq[0]*rdxUx[0]-1104.0*rdxLxCu[0])*rho[2]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[1]+1080.799703922979*rho[0]*rdxUyCu[1]+(71332.78045891662*rho[0]*rdxLy[1]+(66289.04850727606*rdxUx[0]+66289.04850727606*rdxLx[0])*rho[0])*rdxUySq[1]+((65208.24880335309*rdxUxSq[0]+60164.51685171252*rdxLx[0]*rdxUx[0]+65208.24880335309*rdxLxSq[0])*rho[0]-71332.78045891662*rho[0]*rdxLySq[1])*rdxUy[1]-1080.799703922979*rho[0]*rdxLyCu[1]+((-66289.04850727606*rdxUx[0])-66289.04850727606*rdxLx[0])*rho[0]*rdxLySq[1]+((-65208.24880335309*rdxUxSq[0])-60164.51685171252*rdxLx[0]*rdxUx[0]-65208.24880335309*rdxLxSq[0])*rho[0]*rdxLy[1])*omega*volFac+(((9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxUyCu[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLy[1]+9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1]-181.8653347947321*rdxUxCu[0]-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(233.8268590217983*rdxUx[0]*rdxUyCu[1]+(15900.22641348229*rdxUx[0]*rdxLy[1]+10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15900.22641348229*rdxUx[0]*rdxLySq[1]+(47804.60228890101*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0])*rdxLy[1]+10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+415.6921938165305*rdxUxR4[0]+28630.79984911354*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+((-181.8653347947321*rdxUxCu[0])-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15900.22641348229*rdxLx[0]*rdxLy[1])-18758.11024597094*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxUySq[1]+((-15900.22641348229*rdxLx[0]*rdxLySq[1])+((-40893.71956670118*rdxLx[0]*rdxUx[0])-47804.60228890101*rdxLxSq[0])*rdxLy[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxUy[1]-233.8268590217983*rdxLx[0]*rdxLyCu[1]+((-18758.11024597094*rdxLx[0]*rdxUx[0])-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-19303.70625035514*rdxLx[0]*rdxUxSq[0])-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxLy[1]-779.4228634059946*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-28630.79984911354*rdxLxCu[0]*rdxUx[0]-415.6921938165305*rdxLxR4[0])*phiLx[3]+((-396.0*rdxUyR4[1])+((-25758.0*rdxLy[1])-24099.0*rdxUx[0]-24099.0*rdxLx[0])*rdxUyCu[1]+(51462.0*rdxLySq[1]+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLy[1]-23220.0*rdxUxSq[0]+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxUySq[1]+(774.0*rdxLyCu[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1]+483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-243.0*rdxUx[0]*rdxUyCu[1])+((-16524.0*rdxUx[0]*rdxLy[1])-9972.0*rdxUxSq[0]-20322.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-16524.0*rdxUx[0]*rdxLySq[1])+((-50364.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0])*rdxLy[1]-10125.0*rdxUxCu[0]-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-243.0*rdxUx[0]*rdxLyCu[1]+((-9972.0*rdxUxSq[0])-20322.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-10125.0*rdxUxCu[0])-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-396.0*rdxUxR4[0]-27378.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-846.0*rdxLxCu[0]*rdxUx[0])*phiUx[2]+(774.0*rdxLy[1]*rdxUyCu[1]+(51462.0*rdxLySq[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLy[1])*rdxUySq[1]+((-25758.0*rdxLyCu[1])+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-396.0*rdxLyR4[1]+((-24099.0*rdxUx[0])-24099.0*rdxLx[0])*rdxLyCu[1]+((-23220.0*rdxUxSq[0])+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxLySq[1]+(483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-243.0*rdxLx[0]*rdxUyCu[1])+((-16524.0*rdxLx[0]*rdxLy[1])-20322.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*rdxUySq[1]+((-16524.0*rdxLx[0]*rdxLySq[1])+((-41814.0*rdxLx[0]*rdxUx[0])-50364.0*rdxLxSq[0])*rdxLy[1]-20925.0*rdxLx[0]*rdxUxSq[0]-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxUy[1]-243.0*rdxLx[0]*rdxLyCu[1]+((-20322.0*rdxLx[0]*rdxUx[0])-9972.0*rdxLxSq[0])*rdxLySq[1]+((-20925.0*rdxLx[0]*rdxUxSq[0])-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxLy[1]-846.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-27378.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[2]+(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0])*phiC[2]+374.1229744348773*phiUy[0]*rdxUyR4[1]+((24224.46259465831*phiUy[0]+841.7766924784738*phiLy[0])*rdxLy[1]+(8865.0*rdxLx[0]-8865.0*rdxUx[0])*phiUy[1]-585.0*rdxUx[0]*phiUx[1]+585.0*rdxLx[0]*phiLx[1]+(22712.38223965068*phiUy[0]+607.9498334566757*phiUx[0])*rdxUx[0]+(22712.38223965068*phiUy[0]+607.9498334566757*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((56024.91542162288*phiLy[0]-56024.91542162288*phiUy[0])*rdxLySq[1]+((4275.0*rdxLx[0]-4275.0*rdxUx[0])*phiUy[1]-38610.0*rdxUx[0]*phiUx[1]+(19215.0*rdxLx[0]-19215.0*rdxUx[0])*phiLy[1]+38610.0*rdxLx[0]*phiLx[1]+((-44349.16092780109*phiUy[0])+40124.68900814059*phiUx[0]+51862.79733103487*phiLy[0])*rdxUx[0]+((-44349.16092780109*phiUy[0])+51862.79733103487*phiLy[0]+40124.68900814059*phiLx[0])*rdxLx[0])*rdxLy[1]+(8640.0*rdxLxSq[0]-8640.0*rdxUxSq[0])*phiUy[1]+((-24960.0*rdxUxSq[0])-46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(46800.0*rdxLx[0]*rdxUx[0]+24960.0*rdxLxSq[0])*phiLx[1]+(21740.70173660454*phiUy[0]+23777.59348630554*phiUx[0])*rdxUxSq[0]+((-34876.57506120691*phiUy[0])+50797.58608438003*phiUx[0]+50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+(21740.70173660454*phiUy[0]+23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-841.7766924784738*phiUy[0])-24224.46259465831*phiLy[0])*rdxLyCu[1]+((19215.0*rdxUx[0]-19215.0*rdxLx[0])*phiUy[1]+38610.0*rdxUx[0]*phiUx[1]+(4275.0*rdxUx[0]-4275.0*rdxLx[0])*phiLy[1]-38610.0*rdxLx[0]*phiLx[1]+((-51862.79733103487*phiUy[0])-40124.68900814059*phiUx[0]+44349.16092780109*phiLy[0])*rdxUx[0]+((-51862.79733103487*phiUy[0])+44349.16092780109*phiLy[0]-40124.68900814059*phiLx[0])*rdxLx[0])*rdxLySq[1]+((19440.0*rdxUxSq[0]-19440.0*rdxLxSq[0])*phiUy[1]+(19440.0*rdxLxSq[0]-19440.0*rdxUxSq[0])*phiLy[1]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxUxSq[0]+(102561.6565193835*phiLy[0]-102561.6565193835*phiUy[0])*rdxLx[0]*rdxUx[0]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxLxSq[0])*rdxLy[1]+(225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*phiUy[1]+((-24375.0*rdxUxCu[0])-32565.0*rdxLx[0]*rdxUxSq[0]-46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(46215.0*rdxLx[0]*rdxUxSq[0]+32565.0*rdxLxSq[0]*rdxUx[0]+24375.0*rdxLxCu[0])*phiLx[1]+(23169.64365284887*phiUx[0]-597.5575286112626*phiUy[0])*rdxUxCu[0]+((-40633.91194556586*phiUy[0])+33842.54072908828*phiUx[0]+50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40633.91194556586*phiUy[0])+50189.63625092335*phiUx[0]+33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(23169.64365284887*phiLx[0]-597.5575286112626*phiUy[0])*rdxLxCu[0])*rdxUy[1]-374.1229744348773*phiLy[0]*rdxLyR4[1]+(585.0*rdxUx[0]*phiUx[1]+(8865.0*rdxUx[0]-8865.0*rdxLx[0])*phiLy[1]-585.0*rdxLx[0]*phiLx[1]+((-607.9498334566757*phiUx[0])-22712.38223965068*phiLy[0])*rdxUx[0]+((-22712.38223965068*phiLy[0])-607.9498334566757*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((24960.0*rdxUxSq[0]+46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(8640.0*rdxUxSq[0]-8640.0*rdxLxSq[0])*phiLy[1]+((-46800.0*rdxLx[0]*rdxUx[0])-24960.0*rdxLxSq[0])*phiLx[1]+((-23777.59348630554*phiUx[0])-21740.70173660454*phiLy[0])*rdxUxSq[0]+((-50797.58608438003*phiUx[0])+34876.57506120691*phiLy[0]-50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+((-21740.70173660454*phiLy[0])-23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((24375.0*rdxUxCu[0]+32565.0*rdxLx[0]*rdxUxSq[0]+46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*phiLy[1]+((-46215.0*rdxLx[0]*rdxUxSq[0])-32565.0*rdxLxSq[0]*rdxUx[0]-24375.0*rdxLxCu[0])*phiLx[1]+(597.5575286112626*phiLy[0]-23169.64365284887*phiUx[0])*rdxUxCu[0]+((-33842.54072908828*phiUx[0])+40633.91194556586*phiLy[0]-50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-50189.63625092335*phiUx[0])+40633.91194556586*phiLy[0]-33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(597.5575286112626*phiLy[0]-23169.64365284887*phiLx[0])*rdxLxCu[0])*rdxLy[1])*omega+((-36.0*rdxUyR4[1])+((-4824.0*rdxLy[1])-3072.0*rdxUx[0]-3072.0*rdxLx[0])*rdxUyCu[1]+((-161676.0*rdxLySq[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*rdxLy[1]-6072.0*rdxUxSq[0]-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*rdxUySq[1]+((-4824.0*rdxLyCu[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*rdxLySq[1]+((-207144.0*rdxUxSq[0])-500088.0*rdxLx[0]*rdxUx[0]-207144.0*rdxLxSq[0])*rdxLy[1]-3072.0*rdxUxCu[0]-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*rdxUy[1]-36.0*rdxLyR4[1]+((-3072.0*rdxUx[0])-3072.0*rdxLx[0])*rdxLyCu[1]+((-6072.0*rdxUxSq[0])-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*rdxLySq[1]+((-3072.0*rdxUxCu[0])-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*rdxLy[1]-36.0*rdxUxR4[0]-4824.0*rdxLx[0]*rdxUxCu[0]-161676.0*rdxLxSq[0]*rdxUxSq[0]-4824.0*rdxLxCu[0]*rdxUx[0]-36.0*rdxLxR4[0])*phiC[2]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[3] = (((144.0*rdxUyCu[1]+(9792.0*rdxLy[1]+3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxUySq[1]+(9792.0*rdxLySq[1]+(31568.0*rdxUx[0]+31568.0*rdxLx[0])*rdxLy[1]+3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxUy[1]+144.0*rdxLyCu[1]+(3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxLySq[1]+(3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxLy[1]+144.0*rdxUxCu[0]+9792.0*rdxLx[0]*rdxUxSq[0]+9792.0*rdxLxSq[0]*rdxUx[0]+144.0*rdxLxCu[0])*rho[3]+((8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxUySq[1]+((6845.064791512203*rdxUx[0]-6845.064791512203*rdxLx[0])*rdxLy[1]-8646.397631383834*rdxUxSq[0]+8646.397631383834*rdxLxSq[0])*rdxUy[1]+(8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxLySq[1]+(8646.397631383834*rdxLxSq[0]-8646.397631383834*rdxUxSq[0])*rdxLy[1]-360.2665679743264*rdxUxCu[0]-23777.59348630554*rdxLx[0]*rdxUxSq[0]+23777.59348630554*rdxLxSq[0]*rdxUx[0]+360.2665679743264*rdxLxCu[0])*rho[2]+((-360.2665679743264*rdxUyCu[1])+((-23777.59348630554*rdxLy[1])-8646.397631383834*rdxUx[0]-8646.397631383834*rdxLx[0])*rdxUySq[1]+(23777.59348630554*rdxLySq[1]-8286.131063409508*rdxUxSq[0]+6845.064791512203*rdxLx[0]*rdxUx[0]-8286.131063409508*rdxLxSq[0])*rdxUy[1]+360.2665679743264*rdxLyCu[1]+(8646.397631383834*rdxUx[0]+8646.397631383834*rdxLx[0])*rdxLySq[1]+(8286.131063409508*rdxUxSq[0]-6845.064791512203*rdxLx[0]*rdxUx[0]+8286.131063409508*rdxLxSq[0])*rdxLy[1])*rho[1]+(21632.0*rdxUx[0]-21632.0*rdxLx[0])*rho[0]*rdxUySq[1]+(21632.0*rdxUxSq[0]-21632.0*rdxLxSq[0])*rho[0]*rdxUy[1]+(21632.0*rdxLx[0]-21632.0*rdxUx[0])*rho[0]*rdxLySq[1]+(21632.0*rdxLxSq[0]-21632.0*rdxUxSq[0])*rho[0]*rdxLy[1])*omega*volFac+((132.0*rdxUyR4[1]+(8586.0*rdxLy[1]+3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxUyCu[1]+((-17154.0*rdxLySq[1])+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLy[1]+2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-258.0*rdxLyCu[1])+((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1]-63.0*rdxUxCu[0]-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-63.0*rdxUx[0]*rdxUyCu[1])+((-4284.0*rdxUx[0]*rdxLy[1])+2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-4284.0*rdxUx[0]*rdxLySq[1])+((-17516.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0])*rdxLy[1]+3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-63.0*rdxUx[0]*rdxLyCu[1]+(2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+132.0*rdxUxR4[0]+8586.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]-258.0*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((-258.0*rdxLy[1]*rdxUyCu[1])+(((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLy[1]-17154.0*rdxLySq[1])*rdxUySq[1]+(8586.0*rdxLyCu[1]+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+132.0*rdxLyR4[1]+(3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxLyCu[1]+(2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxLySq[1]+((-63.0*rdxUxCu[0])-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-63.0*rdxLx[0]*rdxUyCu[1])+((-4284.0*rdxLx[0]*rdxLy[1])-6158.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-4284.0*rdxLx[0]*rdxLySq[1])+((-10106.0*rdxLx[0]*rdxUx[0])-17516.0*rdxLxSq[0])*rdxLy[1]-6353.0*rdxLx[0]*rdxUxSq[0]-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxUy[1]-63.0*rdxLx[0]*rdxLyCu[1]+(2812.0*rdxLxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6353.0*rdxLx[0]*rdxUxSq[0])-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxLy[1]-258.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]+8586.0*rdxLxCu[0]*rdxUx[0]+132.0*rdxLxR4[0])*phiLx[3]+((-12.0*rdxUyR4[1])+((-1608.0*rdxLy[1])-1024.0*rdxUx[0]-1024.0*rdxLx[0])*rdxUyCu[1]+((-53892.0*rdxLySq[1])+((-69632.0*rdxUx[0])-69632.0*rdxLx[0])*rdxLy[1]-2024.0*rdxUxSq[0]-69048.0*rdxLx[0]*rdxUx[0]-2024.0*rdxLxSq[0])*rdxUySq[1]+((-1608.0*rdxLyCu[1])+((-69632.0*rdxUx[0])-69632.0*rdxLx[0])*rdxLySq[1]+((-69048.0*rdxUxSq[0])-166696.0*rdxLx[0]*rdxUx[0]-69048.0*rdxLxSq[0])*rdxLy[1]-1024.0*rdxUxCu[0]-69632.0*rdxLx[0]*rdxUxSq[0]-69632.0*rdxLxSq[0]*rdxUx[0]-1024.0*rdxLxCu[0])*rdxUy[1]-12.0*rdxLyR4[1]+((-1024.0*rdxUx[0])-1024.0*rdxLx[0])*rdxLyCu[1]+((-2024.0*rdxUxSq[0])-69048.0*rdxLx[0]*rdxUx[0]-2024.0*rdxLxSq[0])*rdxLySq[1]+((-1024.0*rdxUxCu[0])-69632.0*rdxLx[0]*rdxUxSq[0]-69632.0*rdxLxSq[0]*rdxUx[0]-1024.0*rdxLxCu[0])*rdxLy[1]-12.0*rdxUxR4[0]-1608.0*rdxLx[0]*rdxUxCu[0]-53892.0*rdxLxSq[0]*rdxUxSq[0]-1608.0*rdxLxCu[0]*rdxUx[0]-12.0*rdxLxR4[0])*phiC[3]+((8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxUyCu[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLy[1]-7925.864495435182*rdxUxSq[0]+7925.864495435182*rdxLxSq[0])*rdxUySq[1]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1]+157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(77.94228634059945*rdxUx[0]*rdxUyCu[1]+(5300.075471160763*rdxUx[0]*rdxLy[1]-2591.14800812304*rdxUxSq[0]+6730.749438212657*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(5300.075471160763*rdxUx[0]*rdxLySq[1]+(20937.03016189259*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0])*rdxLy[1]-2793.797952608599*rdxUxCu[0]+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+77.94228634059945*rdxUx[0]*rdxLyCu[1]+(6730.749438212657*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxUxSq[0])*rdxLySq[1]+((-2793.797952608599*rdxUxCu[0])+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-124.7076581449591*rdxUxR4[0]-8074.820864886104*rdxLx[0]*rdxUxCu[0]+18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+280.592230826158*rdxLxCu[0]*rdxUx[0])*phiUx[2]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxLyCu[1]+(7925.864495435182*rdxLxSq[0]-7925.864495435182*rdxUxSq[0])*rdxLySq[1]+(157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-77.94228634059945*rdxLx[0]*rdxUyCu[1])+((-5300.075471160763*rdxLx[0]*rdxLy[1])-6730.749438212657*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*rdxUySq[1]+((-5300.075471160763*rdxLx[0]*rdxLySq[1])+((-13236.33227144136*rdxLx[0]*rdxUx[0])-20937.03016189259*rdxLxSq[0])*rdxLy[1]-6933.399382698215*rdxLx[0]*rdxUxSq[0]-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxUy[1]-77.94228634059945*rdxLx[0]*rdxLyCu[1]+(2591.14800812304*rdxLxSq[0]-6730.749438212657*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6933.399382698215*rdxLx[0]*rdxUxSq[0])-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxLy[1]-280.592230826158*rdxLx[0]*rdxUxCu[0]-18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+8074.820864886104*rdxLxCu[0]*rdxUx[0]+124.7076581449591*rdxLxR4[0])*phiLx[2]-124.7076581449591*phiUy[1]*rdxUyR4[1]+(((-8074.820864886104*phiUy[1])-280.592230826158*phiLy[1])*rdxLy[1]+((-2793.797952608599*rdxUx[0])-2793.797952608599*rdxLx[0])*phiUy[1]+157.6166234887678*rdxUx[0]*phiUx[1]+157.6166234887678*rdxLx[0]*phiLx[1]+(7683.0*phiUy[0]-195.0*phiUx[0])*rdxUx[0]+(195.0*phiLx[0]-7683.0*phiUy[0])*rdxLx[0])*rdxUyCu[1]+((18674.97180720763*phiUy[1]-18674.97180720763*phiLy[1])*rdxLySq[1]+((17086.68121666697*rdxUx[0]+17086.68121666697*rdxLx[0])*phiUy[1]+10402.69715025867*rdxUx[0]*phiUx[1]+((-6933.399382698215*rdxUx[0])-6933.399382698215*rdxLx[0])*phiLy[1]+10402.69715025867*rdxLx[0]*phiLx[1]+(3705.0*phiUy[0]-12870.0*phiUx[0]+16653.0*phiLy[0])*rdxUx[0]+((-3705.0*phiUy[0])-16653.0*phiLy[0]+12870.0*phiLx[0])*rdxLx[0])*rdxLy[1]+((-2591.14800812304*rdxUxSq[0])+20937.03016189259*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxLxSq[0])*phiUy[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxUxSq[0])*phiUx[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxLxSq[0])*phiLx[1]+(7488.0*phiUy[0]+7488.0*phiUx[0])*rdxUxSq[0]+(16848.0*phiLx[0]-16848.0*phiUx[0])*rdxLx[0]*rdxUx[0]+((-7488.0*phiUy[0])-7488.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((280.592230826158*phiUy[1]+8074.820864886104*phiLy[1])*rdxLyCu[1]+((6933.399382698215*rdxUx[0]+6933.399382698215*rdxLx[0])*phiUy[1]-10402.69715025867*rdxUx[0]*phiUx[1]+((-17086.68121666697*rdxUx[0])-17086.68121666697*rdxLx[0])*phiLy[1]-10402.69715025867*rdxLx[0]*phiLx[1]+((-16653.0*phiUy[0])+12870.0*phiUx[0]-3705.0*phiLy[0])*rdxUx[0]+(16653.0*phiUy[0]+3705.0*phiLy[0]-12870.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((6730.749438212657*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0]+6730.749438212657*rdxLxSq[0])*phiUy[1]+((-6730.749438212657*rdxUxSq[0])-13236.33227144136*rdxLx[0]*rdxUx[0]-6730.749438212657*rdxLxSq[0])*phiLy[1]+(16848.0*phiLy[0]-16848.0*phiUy[0])*rdxUxSq[0]+(16848.0*phiUy[0]-16848.0*phiLy[0])*rdxLxSq[0])*rdxLy[1]+(77.94228634059945*rdxUxCu[0]+5300.075471160763*rdxLx[0]*rdxUxSq[0]+5300.075471160763*rdxLxSq[0]*rdxUx[0]+77.94228634059945*rdxLxCu[0])*phiUy[1]+((-8083.48111892395*rdxUxCu[0])-2994.715846286589*rdxLx[0]*rdxUxSq[0]+15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(15333.84579940727*rdxLx[0]*rdxUxSq[0]-2994.715846286589*rdxLxSq[0]*rdxUx[0]-8083.48111892395*rdxLxCu[0])*phiLx[1]+(7683.0*phiUx[0]-195.0*phiUy[0])*rdxUxCu[0]+((-12870.0*phiUy[0])+3705.0*phiUx[0]+16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(12870.0*phiUy[0]-16653.0*phiUx[0]-3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(195.0*phiUy[0]-7683.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+124.7076581449591*phiLy[1]*rdxLyR4[1]+((-157.6166234887678*rdxUx[0]*phiUx[1])+(2793.797952608599*rdxUx[0]+2793.797952608599*rdxLx[0])*phiLy[1]-157.6166234887678*rdxLx[0]*phiLx[1]+(195.0*phiUx[0]-7683.0*phiLy[0])*rdxUx[0]+(7683.0*phiLy[0]-195.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((7925.864495435182*rdxUxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiUx[1]+(2591.14800812304*rdxUxSq[0]-20937.03016189259*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*phiLy[1]+(7925.864495435182*rdxLxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiLx[1]+((-7488.0*phiUx[0])-7488.0*phiLy[0])*rdxUxSq[0]+(16848.0*phiUx[0]-16848.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(7488.0*phiLy[0]+7488.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((8083.48111892395*rdxUxCu[0]+2994.715846286589*rdxLx[0]*rdxUxSq[0]-15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-77.94228634059945*rdxUxCu[0])-5300.075471160763*rdxLx[0]*rdxUxSq[0]-5300.075471160763*rdxLxSq[0]*rdxUx[0]-77.94228634059945*rdxLxCu[0])*phiLy[1]+((-15333.84579940727*rdxLx[0]*rdxUxSq[0])+2994.715846286589*rdxLxSq[0]*rdxUx[0]+8083.48111892395*rdxLxCu[0])*phiLx[1]+(195.0*phiLy[0]-7683.0*phiUx[0])*rdxUxCu[0]+((-3705.0*phiUx[0])+12870.0*phiLy[0]-16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(16653.0*phiUx[0]-12870.0*phiLy[0]+3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(7683.0*phiLx[0]-195.0*phiLy[0])*rdxLxCu[0])*rdxLy[1])*omega+(12.0*rdxUyR4[1]+(1608.0*rdxLy[1]+1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxUyCu[1]+(53892.0*rdxLySq[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLy[1]+2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxUySq[1]+(1608.0*rdxLyCu[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLySq[1]+(69048.0*rdxUxSq[0]+166696.0*rdxLx[0]*rdxUx[0]+69048.0*rdxLxSq[0])*rdxLy[1]+1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxUy[1]+12.0*rdxLyR4[1]+(1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxLyCu[1]+(2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxLySq[1]+(1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxLy[1]+12.0*rdxUxR4[0]+1608.0*rdxLx[0]*rdxUxCu[0]+53892.0*rdxLxSq[0]*rdxUxSq[0]+1608.0*rdxLxCu[0]*rdxUx[0]+12.0*rdxLxR4[0])*phiC[3])/(12.0*rdxUyR4[1]+(1608.0*rdxLy[1]+1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxUyCu[1]+(53892.0*rdxLySq[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLy[1]+2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxUySq[1]+(1608.0*rdxLyCu[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLySq[1]+(69048.0*rdxUxSq[0]+166696.0*rdxLx[0]*rdxUx[0]+69048.0*rdxLxSq[0])*rdxLy[1]+1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxUy[1]+12.0*rdxLyR4[1]+(1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxLyCu[1]+(2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxLySq[1]+(1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxLy[1]+12.0*rdxUxR4[0]+1608.0*rdxLx[0]*rdxUxCu[0]+53892.0*rdxLxSq[0]*rdxUxSq[0]+1608.0*rdxLxCu[0]*rdxUx[0]+12.0*rdxLxR4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[0]*rho[1]+144.0*rho[0]*rdxCp2[1]+1600.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(81.0*phiUy[0]+81.0*phiLy[0]-162.0*phiC[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+93.53074360871933*rdxCp2[0]*phiUx[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+(900.0*phiUy[0]-54.0*phiUx[0]+900.0*phiLy[0]-2178.0*phiC[0]+864.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiUx[1]+(1020.0*phiUx[0]-2580.0*phiC[0]+3120.0*bcVals[0])*rdxCp2Sq[0])*omega+162.0*phiC[0]*rdxCp2Sq[1]+2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+2580.0*phiC[0]*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]+277.1281292110203*rdxCp2[0]*rho[0])*omega*volFac+(((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1]-162.0*phiC[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]-360.0*rdxCp2[0]*phiUx[1]+189.0*rdxCp2[0]*phiLy[1]-2178.0*rdxCp2[0]*phiC[1]+(155.8845726811989*phiUy[0]+311.7691453623978*phiUx[0]+155.8845726811989*phiLy[0]-1247.076581449591*bcVals[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiUx[1]-2580.0*rdxCp2Sq[0]*phiC[1]+(623.5382907247956*phiUx[0]-1247.076581449591*bcVals[0])*rdxCp2Sq[0])*omega+162.0*phiC[1]*rdxCp2Sq[1]+2178.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+2580.0*rdxCp2Sq[0]*phiC[1])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = ((748.2459488697547*rdxCp2[0]*rho[3]+(368.0*rdxCp2[1]+1600.0*rdxCp2[0])*rho[2])*omega*volFac+((-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3])+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(1020.0*rdxCp2Sq[0]-138.0*rdxCp2[0]*rdxCp2[1])*phiUx[2]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[2]+(199.1858428704209*phiUy[0]-199.1858428704209*phiLy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiUy[0]-866.0254037844386*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[2])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[0]*rho[2])*omega*volFac+(((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiUx[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+623.5382907247956*rdxCp2Sq[0])*phiUx[2]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiUy[0]-150.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[3])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[1]-576.0*rho[0]*rdxCp2[1]-1680.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-324.0*phiUy[0])-324.0*phiLy[0]+648.0*phiC[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+311.7691453623978*rdxCp2[0]*phiUx[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+((-945.0*phiUy[0])-324.0*phiUx[0]-945.0*phiLy[0]+2214.0*phiC[0]+576.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiUx[1]+((-720.0*phiUx[0])+720.0*phiC[0]+2080.0*bcVals[0])*rdxCp2Sq[0])*omega-648.0*phiC[0]*rdxCp2Sq[1]-2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-720.0*phiC[0]*rdxCp2Sq[0]))/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]-138.5640646055102*rdxCp2[0]*rho[0])*omega*volFac+(((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]+75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1]-216.0*phiC[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiUx[1]+54.0*rdxCp2[0]*phiLy[1]-738.0*rdxCp2[0]*phiC[1]+((-77.94228634059945*phiUy[0])+155.8845726811989*phiUx[0]-77.94228634059945*phiLy[0]+277.1281292110203*bcVals[0])*rdxCp2[0])*rdxCp2[1]-240.0*rdxCp2Sq[0]*phiC[1]+277.1281292110203*bcVals[0]*rdxCp2Sq[0])*omega+216.0*phiC[1]*rdxCp2Sq[1]+738.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+240.0*rdxCp2Sq[0]*phiC[1])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[3]+((-1472.0*rdxCp2[1])-1680.0*rdxCp2[0])*rho[2])*omega*volFac+((-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3])+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-828.0*rdxCp2[0]*rdxCp2[1])-720.0*rdxCp2Sq[0])*phiUx[2]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[2]+(796.7433714816835*phiLy[0]-796.7433714816835*phiUy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiLy[0]-909.3266739736605*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[2]))/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[0]*rho[2])*omega*volFac+(((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[3]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]+1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiUx[2]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiLy[0]-225.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[3])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[1]-144.0*rho[0]*rdxCp2[1]-1600.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-81.0*phiUy[0])-81.0*phiLy[0]+162.0*phiC[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+93.53074360871933*rdxCp2[0]*phiLx[1]-864.0*rdxCp2[0]*bcVals[1]+((-900.0*phiUy[0])-900.0*phiLy[0]+54.0*phiLx[0]+2178.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiLx[1]-3120.0*rdxCp2Sq[0]*bcVals[1]+(2580.0*phiC[0]-1020.0*phiLx[0])*rdxCp2Sq[0])*omega-162.0*phiC[0]*rdxCp2Sq[1]-2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-2580.0*phiC[0]*rdxCp2Sq[0]))/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]-277.1281292110203*rdxCp2[0]*rho[0])*omega*volFac+(((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1]-162.0*phiC[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]+189.0*rdxCp2[0]*phiLy[1]-360.0*rdxCp2[0]*phiLx[1]-2178.0*rdxCp2[0]*phiC[1]+1247.076581449591*rdxCp2[0]*bcVals[1]+((-155.8845726811989*phiUy[0])-155.8845726811989*phiLy[0]-311.7691453623978*phiLx[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiLx[1]-2580.0*rdxCp2Sq[0]*phiC[1]+1247.076581449591*rdxCp2Sq[0]*bcVals[1]-623.5382907247956*phiLx[0]*rdxCp2Sq[0])*omega+162.0*phiC[1]*rdxCp2Sq[1]+2178.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+2580.0*rdxCp2Sq[0]*phiC[1])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[3]+((-368.0*rdxCp2[1])-1600.0*rdxCp2[0])*rho[2])*omega*volFac+((-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3])-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiLx[3]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(138.0*rdxCp2[0]*rdxCp2[1]-1020.0*rdxCp2Sq[0])*phiLx[2]+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[2]+(199.1858428704209*phiLy[0]-199.1858428704209*phiUy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiLy[0]-866.0254037844386*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[2]))/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[0]*rho[2])*omega*volFac+(((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiLx[3]+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[3]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-796.7433714816835*rdxCp2[0]*rdxCp2[1])-623.5382907247956*rdxCp2Sq[0])*phiLx[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiLy[0]-150.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[3])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[0]*rho[1]+576.0*rho[0]*rdxCp2[1]+1680.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(324.0*phiUy[0]+324.0*phiLy[0]-648.0*phiC[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+311.7691453623978*rdxCp2[0]*phiLx[1]+576.0*rdxCp2[0]*bcVals[1]+(945.0*phiUy[0]+945.0*phiLy[0]+324.0*phiLx[0]-2214.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiLx[1]+2080.0*rdxCp2Sq[0]*bcVals[1]+(720.0*phiLx[0]-720.0*phiC[0])*rdxCp2Sq[0])*omega+648.0*phiC[0]*rdxCp2Sq[1]+2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+720.0*phiC[0]*rdxCp2Sq[0])/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]+138.5640646055102*rdxCp2[0]*rho[0])*omega*volFac+(((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]-75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1]-216.0*phiC[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]+54.0*rdxCp2[0]*phiLy[1]-150.0*rdxCp2[0]*phiLx[1]-738.0*rdxCp2[0]*phiC[1]+277.1281292110203*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiUy[0]+77.94228634059945*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1]-240.0*rdxCp2Sq[0]*phiC[1]+277.1281292110203*rdxCp2Sq[0]*bcVals[1])*omega+216.0*phiC[1]*rdxCp2Sq[1]+738.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+240.0*rdxCp2Sq[0]*phiC[1])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = ((277.1281292110203*rdxCp2[0]*rho[3]+(1472.0*rdxCp2[1]+1680.0*rdxCp2[0])*rho[2])*omega*volFac+((-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3])-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiLx[3]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(828.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiLx[2]+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[2]+(796.7433714816835*phiUy[0]-796.7433714816835*phiLy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiUy[0]-909.3266739736605*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[2])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[0]*rho[2])*omega*volFac+(((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[3]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]-1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiUy[0]-225.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[3])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[1]*rho[2]+1600.0*rho[0]*rdxCp2[1]+144.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(3120.0*rdxCp2Sq[1]+864.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(1020.0*phiUy[0]-2580.0*phiC[0])*rdxCp2Sq[1]+((-866.0254037844386*rdxCp2[0]*phiUx[1])+866.0254037844386*rdxCp2[0]*phiLx[1]+((-54.0*phiUy[0])+900.0*phiUx[0]+900.0*phiLx[0]-2178.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-77.94228634059945*rdxCp2Sq[0]*phiUx[1]+77.94228634059945*rdxCp2Sq[0]*phiLx[1]+(81.0*phiUx[0]+81.0*phiLx[0]-162.0*phiC[0])*rdxCp2Sq[0])*omega+2580.0*phiC[0]*rdxCp2Sq[1]+2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+162.0*phiC[0]*rdxCp2Sq[0])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = ((748.2459488697547*rdxCp2[1]*rho[3]+(1600.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[1])*omega*volFac+((239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(1020.0*phiUy[1]-2580.0*phiC[1])*rdxCp2Sq[1]+((-138.0*rdxCp2[0]*phiUy[1])-700.0*rdxCp2[0]*phiUx[1]-700.0*rdxCp2[0]*phiLx[1]-5566.0*rdxCp2[0]*phiC[1]+(866.0254037844386*phiUx[0]-866.0254037844386*phiLx[0])*rdxCp2[0])*rdxCp2[1]-161.0*rdxCp2Sq[0]*phiUx[1]-161.0*rdxCp2Sq[0]*phiLx[1]-1058.0*rdxCp2Sq[0]*phiC[1]+(199.1858428704209*phiUx[0]-199.1858428704209*phiLx[0])*rdxCp2Sq[0])*omega+2580.0*phiC[1]*rdxCp2Sq[1]+5566.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+1058.0*rdxCp2Sq[0]*phiC[1])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]+277.1281292110203*rho[0]*rdxCp2[1])*omega*volFac+(((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]+((-2580.0*rdxCp2Sq[1])-2178.0*rdxCp2[0]*rdxCp2[1]-162.0*rdxCp2Sq[0])*phiC[2]+((-1247.076581449591*rdxCp2Sq[1])-1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[2]+623.5382907247956*phiUy[0]*rdxCp2Sq[1]+((-150.0*rdxCp2[0]*phiUx[1])+150.0*rdxCp2[0]*phiLx[1]+(311.7691453623978*phiUy[0]+155.8845726811989*phiUx[0]+155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0])*phiC[2])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[1]*rho[1])*omega*volFac+(((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+((-2580.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-1058.0*rdxCp2Sq[0])*phiC[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]+623.5382907247956*phiUy[1]*rdxCp2Sq[1]+(796.7433714816835*rdxCp2[0]*phiUy[1]-121.2435565298214*rdxCp2[0]*phiUx[1]-121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiUx[0]-150.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0])*phiC[3])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[2]-1680.0*rho[0]*rdxCp2[1]-576.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(720.0*phiC[0]-720.0*phiUy[0])*rdxCp2Sq[1]+(909.3266739736605*rdxCp2[0]*phiUx[1]-909.3266739736605*rdxCp2[0]*phiLx[1]+((-324.0*phiUy[0])-945.0*phiUx[0]-945.0*phiLx[0]+2214.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+311.7691453623978*rdxCp2Sq[0]*phiUx[1]-311.7691453623978*rdxCp2Sq[0]*phiLx[1]+((-324.0*phiUx[0])-324.0*phiLx[0]+648.0*phiC[0])*rdxCp2Sq[0])*omega-720.0*phiC[0]*rdxCp2Sq[1]-2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-648.0*phiC[0]*rdxCp2Sq[0]))/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[3]+((-1680.0*rdxCp2[1])-1472.0*rdxCp2[0])*rho[1])*omega*volFac+((692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiUy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(720.0*phiC[1]-720.0*phiUy[1])*rdxCp2Sq[1]+((-828.0*rdxCp2[0]*phiUy[1])+735.0*rdxCp2[0]*phiUx[1]+735.0*rdxCp2[0]*phiLx[1]+5658.0*rdxCp2[0]*phiC[1]+(909.3266739736605*phiLx[0]-909.3266739736605*phiUx[0])*rdxCp2[0])*rdxCp2[1]+644.0*rdxCp2Sq[0]*phiUx[1]+644.0*rdxCp2Sq[0]*phiLx[1]+4232.0*rdxCp2Sq[0]*phiC[1]+(796.7433714816835*phiLx[0]-796.7433714816835*phiUx[0])*rdxCp2Sq[0])*omega-720.0*phiC[1]*rdxCp2Sq[1]-5658.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-4232.0*rdxCp2Sq[0]*phiC[1]))/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]-138.5640646055102*rho[0]*rdxCp2[1])*omega*volFac+(((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+((-240.0*rdxCp2Sq[1])-738.0*rdxCp2[0]*rdxCp2[1]-216.0*rdxCp2Sq[0])*phiC[2]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(75.0*rdxCp2[0]*phiUx[1]-75.0*rdxCp2[0]*phiLx[1]+(155.8845726811989*phiUy[0]-77.94228634059945*phiUx[0]-77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0])*phiC[2])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[1]*rho[1])*omega*volFac+((-1150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+((-720.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-4232.0*rdxCp2Sq[0])*phiC[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+(1195.115057222525*rdxCp2[0]*phiUy[1]+181.8653347947321*rdxCp2[0]*phiUx[1]+181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiLx[0]-225.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0])*phiC[3])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[2]-1600.0*rho[0]*rdxCp2[1]-144.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-3120.0*rdxCp2Sq[1])-864.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2580.0*phiC[0]-1020.0*phiLy[0])*rdxCp2Sq[1]+(866.0254037844386*rdxCp2[0]*phiUx[1]-866.0254037844386*rdxCp2[0]*phiLx[1]+((-900.0*phiUx[0])+54.0*phiLy[0]-900.0*phiLx[0]+2178.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0]*phiUx[1]-77.94228634059945*rdxCp2Sq[0]*phiLx[1]+((-81.0*phiUx[0])-81.0*phiLx[0]+162.0*phiC[0])*rdxCp2Sq[0])*omega-2580.0*phiC[0]*rdxCp2Sq[1]-2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-162.0*phiC[0]*rdxCp2Sq[0]))/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[3]+((-1600.0*rdxCp2[1])-368.0*rdxCp2[0])*rho[1])*omega*volFac+((-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3])+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2580.0*phiC[1]-1020.0*phiLy[1])*rdxCp2Sq[1]+(700.0*rdxCp2[0]*phiUx[1]+138.0*rdxCp2[0]*phiLy[1]+700.0*rdxCp2[0]*phiLx[1]+5566.0*rdxCp2[0]*phiC[1]+(866.0254037844386*phiLx[0]-866.0254037844386*phiUx[0])*rdxCp2[0])*rdxCp2[1]+161.0*rdxCp2Sq[0]*phiUx[1]+161.0*rdxCp2Sq[0]*phiLx[1]+1058.0*rdxCp2Sq[0]*phiC[1]+(199.1858428704209*phiLx[0]-199.1858428704209*phiUx[0])*rdxCp2Sq[0])*omega-2580.0*phiC[1]*rdxCp2Sq[1]-5566.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-1058.0*rdxCp2Sq[0]*phiC[1]))/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]-277.1281292110203*rho[0]*rdxCp2[1])*omega*volFac+(((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+(1247.076581449591*rdxCp2Sq[1]+1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]+((-2580.0*rdxCp2Sq[1])-2178.0*rdxCp2[0]*rdxCp2[1]-162.0*rdxCp2Sq[0])*phiC[2]-623.5382907247956*phiLy[0]*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiLx[1]+((-155.8845726811989*phiUx[0])-311.7691453623978*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0])*phiC[2])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[1]*rho[1])*omega*volFac+(((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+((-2580.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-1058.0*rdxCp2Sq[0])*phiC[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]-623.5382907247956*phiLy[1]*rdxCp2Sq[1]+(121.2435565298214*rdxCp2[0]*phiUx[1]-796.7433714816835*rdxCp2[0]*phiLy[1]+121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiLx[0]-150.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0])*phiC[3])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[1]*rho[2]+1680.0*rho[0]*rdxCp2[1]+576.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(720.0*phiLy[0]-720.0*phiC[0])*rdxCp2Sq[1]+((-909.3266739736605*rdxCp2[0]*phiUx[1])+909.3266739736605*rdxCp2[0]*phiLx[1]+(945.0*phiUx[0]+324.0*phiLy[0]+945.0*phiLx[0]-2214.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-311.7691453623978*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLx[1]+(324.0*phiUx[0]+324.0*phiLx[0]-648.0*phiC[0])*rdxCp2Sq[0])*omega+720.0*phiC[0]*rdxCp2Sq[1]+2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+648.0*phiC[0]*rdxCp2Sq[0])/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = ((277.1281292110203*rdxCp2[1]*rho[3]+(1680.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[1])*omega*volFac+((-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3])+(692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiLy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(720.0*phiLy[1]-720.0*phiC[1])*rdxCp2Sq[1]+((-735.0*rdxCp2[0]*phiUx[1])+828.0*rdxCp2[0]*phiLy[1]-735.0*rdxCp2[0]*phiLx[1]-5658.0*rdxCp2[0]*phiC[1]+(909.3266739736605*phiUx[0]-909.3266739736605*phiLx[0])*rdxCp2[0])*rdxCp2[1]-644.0*rdxCp2Sq[0]*phiUx[1]-644.0*rdxCp2Sq[0]*phiLx[1]-4232.0*rdxCp2Sq[0]*phiC[1]+(796.7433714816835*phiUx[0]-796.7433714816835*phiLx[0])*rdxCp2Sq[0])*omega+720.0*phiC[1]*rdxCp2Sq[1]+5658.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+4232.0*rdxCp2Sq[0]*phiC[1])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]+138.5640646055102*rho[0]*rdxCp2[1])*omega*volFac+(((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+((-240.0*rdxCp2Sq[1])-738.0*rdxCp2[0]*rdxCp2[1]-216.0*rdxCp2Sq[0])*phiC[2]+((-75.0*rdxCp2[0]*phiUx[1])+75.0*rdxCp2[0]*phiLx[1]+(77.94228634059945*phiUx[0]-155.8845726811989*phiLy[0]+77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0])*phiC[2])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[1]*rho[1])*omega*volFac+(((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+((-720.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-4232.0*rdxCp2Sq[0])*phiC[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+((-181.8653347947321*rdxCp2[0]*phiUx[1])-1195.115057222525*rdxCp2[0]*phiLy[1]-181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiUx[0]-225.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0])*phiC[3])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(328950.0*phiUy[0]-832050.0*phiC[0])*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiUx[1]+(2715915.0*phiUy[0]-193500.0*phiUx[0]-8611395.0*phiC[0]+3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiUx[1]+((-193500.0*phiUy[0])+2715915.0*phiUx[0]-8611395.0*phiC[0]+9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiUx[1]+(328950.0*phiUx[0]-832050.0*phiC[0]+1006200.0*bcVals[0])*rdxCp2R3[0])*omega+832050.0*phiC[0]*rdxCp2R3[1]+8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+832050.0*phiC[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(328950.0*phiUy[1]-832050.0*phiC[1])*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiUx[1]-8611395.0*rdxCp2[0]*phiC[1]+(567939.4598018348*phiUy[0]+1117172.770881926*phiUx[0]-4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiUx[1]-8611395.0*rdxCp2Sq[0]*phiC[1]+((-33515.18312645776*phiUy[0])+1919718.512568964*phiUx[0]-4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiUx[1]-832050.0*rdxCp2R3[0]*phiC[1]+(201091.0987587466*phiUx[0]-402182.1975174932*bcVals[0])*rdxCp2R3[0])*omega+832050.0*phiC[1]*rdxCp2R3[1]+8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+832050.0*rdxCp2R3[0]*phiC[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiUx[1]+(1919718.512568964*phiUy[0]-33515.18312645776*phiUx[0]+536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiUx[1]+(1117172.770881926*phiUy[0]+567939.4598018348*phiUx[0]+1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiUx[1]+(44400.0*phiUy[0]+64500.0*phiUx[0]-258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiUx[1]+(64500.0*phiUy[0]+44400.0*phiUx[0]+40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiC[0]-12960.0*phiUy[0])*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiUx[1]+((-456408.0*phiUy[0])+11340.0*phiUx[0]+535788.0*phiC[0]-181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiUx[1]+((-1842246.0*phiUy[0])-532953.0*phiUx[0]+3688425.0*phiC[0]-2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiUx[1]+((-928800.0*phiUy[0])-2563956.0*phiUx[0]+7679628.0*phiC[0]-8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiUx[1]+((-1052640.0*phiUx[0])+2662560.0*phiC[0]-3219840.0*bcVals[0])*rdxCp2R4[0])*omega-12960.0*phiC[0]*rdxCp2R4[1]-535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-2662560.0*phiC[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+((12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiC[1]-12960.0*phiUy[1])*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiUx[1]+535788.0*rdxCp2[0]*phiC[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiUx[1]+3688425.0*rdxCp2Sq[0]*phiC[1]+((-301792.5327108011*phiUy[0])-659547.6270141526*phiUx[0]+1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiUx[1]+7679628.0*rdxCp2R3[0]*phiC[1]+((-160872.8790069972*phiUy[0])-1745283.675738703*phiUx[0]+3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiUx[1]+2662560.0*rdxCp2R4[0]*phiC[1]+(1286983.032055978*bcVals[0]-643491.516027989*phiUx[0])*rdxCp2R4[0])*omega-12960.0*phiC[1]*rdxCp2R4[1]-535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-2662560.0*rdxCp2R4[0]*phiC[1]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiUx[1]+(98207.28078915528*phiUy[0]+14029.6115413079*phiUx[0]-224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiUx[1]+(2061183.762277152*phiUy[0]-814886.6036909672*phiUx[0]-2492594.31717237*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiUx[1]+(6703036.625291553*phiUy[0]-3407636.758811008*phiUx[0]-6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiUx[1]+(332100.0*bcVals[0]-166050.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiUx[1]+(387000.0*phiUy[0]-266400.0*phiUx[0]-241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(2662560.0*phiC[0]-1052640.0*phiUy[0])*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiUx[1]+((-2563956.0*phiUy[0])-928800.0*phiUx[0]+7679628.0*phiC[0]+1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiUx[1]+((-532953.0*phiUy[0])-1842246.0*phiUx[0]+3688425.0*phiC[0]+5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiUx[1]+(11340.0*phiUy[0]-456408.0*phiUx[0]+535788.0*phiC[0]+1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiUx[1]+((-12960.0*phiUx[0])+12960.0*phiC[0]+37440.0*bcVals[0])*rdxCp2R4[0])*omega-2662560.0*phiC[0]*rdxCp2R4[1]-7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-12960.0*phiC[0]*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(5263200.0*phiUy[1]-1.33128e+7*phiC[1])*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiUx[1]-3.839814e+7*rdxCp2[0]*phiC[1]+((-3407636.758811008*phiUy[0])+6703036.625291553*phiUx[0]+1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiUx[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]+((-814886.6036909672*phiUy[0])+2061183.762277152*phiUx[0]+1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiUx[1]-2678940.0*rdxCp2R3[0]*phiC[1]+(14029.6115413079*phiUy[0]+98207.28078915528*phiUx[0]+2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]+74824.59488697546*bcVals[0]*rdxCp2R4[0])*omega+1.33128e+7*phiC[1]*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiUx[1]+((-1745283.675738703*phiUy[0])-160872.8790069972*phiUx[0]+285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiUx[1]+((-659547.6270141526*phiUy[0])-301792.5327108011*phiUx[0]+871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiUx[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiUx[1]+((-266400.0*phiUy[0])+387000.0*phiUx[0]+688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+(580800.0*bcVals[0]-166050.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiUy[0]-1440.0*phiC[0])*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiUx[1]+(11799.0*phiUy[0]+1890.0*phiUx[0]-13689.0*phiC[0]-3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiUx[1]+(1890.0*phiUy[0]+11799.0*phiUx[0]-13689.0*phiC[0]-33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiUx[1]+(1440.0*phiUx[0]-1440.0*phiC[0]-4160.0*bcVals[0])*rdxCp2R3[0])*omega+1440.0*phiC[0]*rdxCp2R3[1]+13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+1440.0*phiC[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiC[1]-1440.0*phiUy[1])*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiUx[1]+13689.0*rdxCp2[0]*phiC[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiUx[0]-4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiUx[1]+13689.0*rdxCp2Sq[0]*phiC[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiUx[0]-14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1]-1662.768775266122*bcVals[0]*rdxCp2R3[0])*omega-1440.0*phiC[1]*rdxCp2R3[1]-13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiUx[1]+((-467.6537180435967*phiUy[0])+467.6537180435967*phiUx[0]-831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiUx[1]+((-2727.980021920981*phiUy[0])+2727.980021920981*phiUx[0]-8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiUx[1]+(1125.0*phiUy[0]-1125.0*phiUx[0]-2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiUx[1]+((-1125.0*phiUy[0])+1125.0*phiUx[0]-4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(832050.0*phiC[0]-328950.0*phiLy[0])*rdxCp2R3[1]+((-335151.8312645776*rdxCp2[0]*phiUx[1])-1533436.541464953*rdxCp2[0]*phiLy[1]+(193500.0*phiUx[0]-2715915.0*phiLy[0]+8611395.0*phiC[0]-3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2176434.423012786*rdxCp2Sq[0]*phiUx[1]+90490.99444143593*rdxCp2Sq[0]*phiLy[1]+((-2715915.0*phiUx[0])+193500.0*phiLy[0]+8611395.0*phiC[0]-9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiUx[1]+((-328950.0*phiUx[0])+832050.0*phiC[0]-1006200.0*bcVals[0])*rdxCp2R3[0])*omega-832050.0*phiC[0]*rdxCp2R3[1]-8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-832050.0*phiC[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+(832050.0*phiC[1]-328950.0*phiLy[1])*rdxCp2R3[1]+(1290000.0*rdxCp2[0]*phiUx[1]-125505.0*rdxCp2[0]*phiLy[1]+8611395.0*rdxCp2[0]*phiC[1]+((-1117172.770881926*phiUx[0])-567939.4598018348*phiLy[0]+4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2054550.0*rdxCp2Sq[0]*phiUx[1]+40635.0*rdxCp2Sq[0]*phiLy[1]+8611395.0*rdxCp2Sq[0]*phiC[1]+((-1919718.512568964*phiUx[0])+33515.18312645776*phiLy[0]+4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiUx[1]+832050.0*rdxCp2R3[0]*phiC[1]+(402182.1975174932*bcVals[0]-201091.0987587466*phiUx[0])*rdxCp2R3[0])*omega-832050.0*phiC[1]*rdxCp2R3[1]-8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-832050.0*rdxCp2R3[0]*phiC[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]-201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-58050.0*rdxCp2[0]*phiUx[1])-359640.0*rdxCp2[0]*phiLy[1]+(33515.18312645776*phiUx[0]-1919718.512568964*phiLy[0]-536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(493650.0*rdxCp2Sq[0]*phiUx[1]-522450.0*rdxCp2Sq[0]*phiLy[1]+((-567939.4598018348*phiUx[0])-1117172.770881926*phiLy[0]-1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+(74478.18472546172*rdxCp2[0]*phiUx[1]-437394.7904353685*rdxCp2[0]*phiLy[1]+((-64500.0*phiUx[0])-44400.0*phiLy[0]+258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(40789.79651824706*rdxCp2Sq[0]*phiUx[1]-78202.0939617348*rdxCp2Sq[0]*phiLy[1]+((-44400.0*phiUx[0])-64500.0*phiLy[0]-40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(12960.0*phiLy[0]-12960.0*phiC[0])*rdxCp2R4[1]+(19641.45615783106*rdxCp2[0]*phiUx[1]+176773.1054204795*rdxCp2[0]*phiLy[1]+((-11340.0*phiUx[0])+456408.0*phiLy[0]-535788.0*phiC[0]+181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-386469.0325912283*rdxCp2Sq[0]*phiUx[1])+814839.8383191626*rdxCp2Sq[0]*phiLy[1]+(532953.0*phiUx[0]+1842246.0*phiLy[0]-3688425.0*phiC[0]+2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2064285.865273508*rdxCp2R3[0]*phiUx[1])+434356.7733188925*rdxCp2R3[0]*phiLy[1]+(2563956.0*phiUx[0]+928800.0*phiLy[0]-7679628.0*phiC[0]+8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiUx[1]+(1052640.0*phiUx[0]-2662560.0*phiC[0]+3219840.0*bcVals[0])*rdxCp2R4[0])*omega+12960.0*phiC[0]*rdxCp2R4[1]+535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+2662560.0*phiC[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(12960.0*phiLy[1]-12960.0*phiC[1])*rdxCp2R4[1]+((-75600.0*rdxCp2[0]*phiUx[1])+157788.0*rdxCp2[0]*phiLy[1]-535788.0*rdxCp2[0]*phiC[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-727155.0*rdxCp2Sq[0]*phiUx[1])+465750.0*rdxCp2Sq[0]*phiLy[1]-3688425.0*rdxCp2Sq[0]*phiC[1]+(659547.6270141526*phiUx[0]+301792.5327108011*phiLy[0]-1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1862820.0*rdxCp2R3[0]*phiUx[1])+195048.0*rdxCp2R3[0]*phiLy[1]-7679628.0*rdxCp2R3[0]*phiC[1]+(1745283.675738703*phiUx[0]+160872.8790069972*phiLy[0]-3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiUx[1]-2662560.0*rdxCp2R4[0]*phiC[1]+(643491.516027989*phiUx[0]-1286983.032055978*bcVals[0])*rdxCp2R4[0])*omega+12960.0*phiC[1]*rdxCp2R4[1]+535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+2662560.0*rdxCp2R4[0]*phiC[1])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]+(24300.0*rdxCp2[0]*phiUx[1]+218700.0*rdxCp2[0]*phiLy[1]+((-14029.6115413079*phiUx[0])-98207.28078915528*phiLy[0]+224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((814886.6036909672*phiUx[0]-2061183.762277152*phiLy[0]+2492594.31717237*bcVals[0])*rdxCp2Sq[0]-664200.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-2961900.0*rdxCp2R3[0]*phiUx[1])-3134700.0*rdxCp2R3[0]*phiLy[1]+(3407636.758811008*phiUx[0]-6703036.625291553*phiLy[0]+6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-31176.91453623978*rdxCp2[0]*phiUx[1])-155884.5726811989*rdxCp2[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-175759.8556980518*rdxCp2Sq[0]*phiUx[1])-687061.2540923841*rdxCp2Sq[0]*phiLy[1]+(166050.0*phiUx[0]-332100.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-244738.7791094823*rdxCp2R3[0]*phiUx[1])-469212.5637704087*rdxCp2R3[0]*phiLy[1]+(266400.0*phiUx[0]-387000.0*phiLy[0]+241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+((417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1052640.0*phiLy[0]-2662560.0*phiC[0])*rdxCp2R4[1]+((-893738.2167055405*rdxCp2[0]*phiUx[1])-454351.5678414679*rdxCp2[0]*phiLy[1]+(928800.0*phiUx[0]+2563956.0*phiLy[0]-7679628.0*phiC[0]-1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-1772702.040022519*rdxCp2Sq[0]*phiUx[1])-108651.5471587956*rdxCp2Sq[0]*phiLy[1]+(1842246.0*phiUx[0]+532953.0*phiLy[0]-3688425.0*phiC[0]-5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-439178.8027671645*rdxCp2R3[0]*phiUx[1])+1870.614872174387*rdxCp2R3[0]*phiLy[1]+(456408.0*phiUx[0]-11340.0*phiLy[0]-535788.0*phiC[0]-1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiUx[1]+(12960.0*phiUx[0]-12960.0*phiC[0]-37440.0*bcVals[0])*rdxCp2R4[0])*omega+2662560.0*phiC[0]*rdxCp2R4[1]+7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+12960.0*phiC[0]*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+((210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1.33128e+7*phiC[1]-5263200.0*phiLy[1])*rdxCp2R4[1]+(6450000.0*rdxCp2[0]*phiUx[1]-1408860.0*rdxCp2[0]*phiLy[1]+3.839814e+7*rdxCp2[0]*phiC[1]+((-6703036.625291553*phiUx[0])+3407636.758811008*phiLy[0]-1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(1983375.0*rdxCp2Sq[0]*phiUx[1]+63990.0*rdxCp2Sq[0]*phiLy[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]+((-2061183.762277152*phiUx[0])+814886.6036909672*phiLy[0]-1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(94500.0*rdxCp2R3[0]*phiUx[1]+9720.0*rdxCp2R3[0]*phiLy[1]+2678940.0*rdxCp2R3[0]*phiC[1]+((-98207.28078915528*phiUx[0])-14029.6115413079*phiLy[0]-2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1]-74824.59488697546*bcVals[0]*rdxCp2R4[0])*omega-1.33128e+7*phiC[1]*rdxCp2R4[1]-3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2]+643491.516027989*phiLy[0]*rdxCp2R4[1]+((-154800.0*rdxCp2[0]*phiUx[1])-106560.0*rdxCp2[0]*phiLy[1]+(160872.8790069972*phiUx[0]+1745283.675738703*phiLy[0]-285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-290400.0*rdxCp2Sq[0]*phiUx[1])-66420.0*rdxCp2Sq[0]*phiLy[1]+(301792.5327108011*phiUx[0]+659547.6270141526*phiLy[0]-871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-63000.0*rdxCp2R3[0]*phiUx[1])-10800.0*rdxCp2R3[0]*phiLy[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+(372390.9236273086*rdxCp2[0]*phiUx[1]-2016730.678300898*rdxCp2[0]*phiLy[1]+((-387000.0*phiUx[0])+266400.0*phiLy[0]-688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((166050.0*phiLy[0]-580800.0*bcVals[0])*rdxCp2Sq[0]-543205.7742697513*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-25980.76211353316*rdxCp2R3[0]*phiUx[1])-18706.14872174387*rdxCp2R3[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(1440.0*phiC[0]-1440.0*phiLy[0])*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiUx[1]+1818.653347947321*rdxCp2[0]*phiLy[1]+((-1890.0*phiUx[0])-11799.0*phiLy[0]+13689.0*phiC[0]+3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(11353.59304361399*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLy[1]+((-11799.0*phiUx[0])-1890.0*phiLy[0]+13689.0*phiC[0]+33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiUx[1]+((-1440.0*phiUx[0])+1440.0*phiC[0]+4160.0*bcVals[0])*rdxCp2R3[0])*omega-1440.0*phiC[0]*rdxCp2R3[1]-13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-1440.0*phiC[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(1440.0*phiLy[1]-1440.0*phiC[1])*rdxCp2R3[1]+((-2625.0*rdxCp2[0]*phiUx[1])+2664.0*rdxCp2[0]*phiLy[1]-13689.0*rdxCp2[0]*phiC[1]+(2727.980021920981*phiUx[0]-2727.980021920981*phiLy[0]+4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiUx[1])+324.0*rdxCp2Sq[0]*phiLy[1]-13689.0*rdxCp2Sq[0]*phiC[1]+(467.6537180435967*phiUx[0]-467.6537180435967*phiLy[0]+14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]+1662.768775266122*bcVals[0]*rdxCp2R3[0])*omega+1440.0*phiC[1]*rdxCp2R3[1]+13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2]+(450.0*rdxCp2[0]*phiUx[1]+450.0*rdxCp2[0]*phiLy[1]+((-467.6537180435967*phiUx[0])+467.6537180435967*phiLy[0]+831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2625.0*rdxCp2Sq[0]*phiUx[1]-450.0*rdxCp2Sq[0]*phiLy[1]+((-2727.980021920981*phiUx[0])+2727.980021920981*phiLy[0]+8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-1082.531754730548*rdxCp2[0]*phiUx[1])-4546.633369868302*rdxCp2[0]*phiLy[1]+(1125.0*phiUx[0]-1125.0*phiLy[0]+2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1082.531754730548*rdxCp2Sq[0]*phiUx[1]-779.4228634059946*rdxCp2Sq[0]*phiLy[1]+((-1125.0*phiUx[0])+1125.0*phiLy[0]+4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(832050.0*phiC[0]-328950.0*phiUy[0])*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiLx[1]-3096000.0*rdxCp2[0]*bcVals[1]+((-2715915.0*phiUy[0])+193500.0*phiLx[0]+8611395.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiLx[1]-9081960.0*rdxCp2Sq[0]*bcVals[1]+(193500.0*phiUy[0]-2715915.0*phiLx[0]+8611395.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiLx[1]-1006200.0*rdxCp2R3[0]*bcVals[1]+(832050.0*phiC[0]-328950.0*phiLx[0])*rdxCp2R3[0])*omega-832050.0*phiC[0]*rdxCp2R3[1]-8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-832050.0*phiC[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(328950.0*phiUy[1]-832050.0*phiC[1])*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiLx[1]-8611395.0*rdxCp2[0]*phiC[1]+4468691.083527703*rdxCp2[0]*bcVals[1]+((-567939.4598018348*phiUy[0])-1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiLx[1]-8611395.0*rdxCp2Sq[0]*phiC[1]+4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(33515.18312645776*phiUy[0]-1919718.512568964*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiLx[1]-832050.0*rdxCp2R3[0]*phiC[1]+402182.1975174932*rdxCp2R3[0]*bcVals[1]-201091.0987587466*phiLx[0]*rdxCp2R3[0])*omega+832050.0*phiC[1]*rdxCp2R3[1]+8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+832050.0*rdxCp2R3[0]*phiC[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiLx[1]-536242.9300233242*rdxCp2[0]*bcVals[1]+(33515.18312645776*phiLx[0]-1919718.512568964*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiLx[1]-1098466.622160182*rdxCp2Sq[0]*bcVals[1]+((-1117172.770881926*phiUy[0])-567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiLx[1]+258000.0*rdxCp2[0]*bcVals[1]+((-44400.0*phiUy[0])-64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiLx[1]-40200.0*rdxCp2Sq[0]*bcVals[1]+((-64500.0*phiUy[0])-44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiUy[0]-12960.0*phiC[0])*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiLx[1]+181440.0*rdxCp2[0]*bcVals[1]+(456408.0*phiUy[0]-11340.0*phiLx[0]-535788.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiLx[1]+2626452.0*rdxCp2Sq[0]*bcVals[1]+(1842246.0*phiUy[0]+532953.0*phiLx[0]-3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiLx[1]+8373744.0*rdxCp2R3[0]*bcVals[1]+(928800.0*phiUy[0]+2563956.0*phiLx[0]-7679628.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiLx[1]+3219840.0*rdxCp2R4[0]*bcVals[1]+(1052640.0*phiLx[0]-2662560.0*phiC[0])*rdxCp2R4[0])*omega+12960.0*phiC[0]*rdxCp2R4[1]+535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+2662560.0*phiC[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+((12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiC[1]-12960.0*phiUy[1])*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiLx[1]+535788.0*rdxCp2[0]*phiC[1]-261886.0821044141*rdxCp2[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiLx[1]+3688425.0*rdxCp2Sq[0]*phiC[1]-1922680.319449907*rdxCp2Sq[0]*bcVals[1]+(301792.5327108011*phiUy[0]+659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiLx[1]+7679628.0*rdxCp2R3[0]*phiC[1]-3812313.1094914*rdxCp2R3[0]*bcVals[1]+(160872.8790069972*phiUy[0]+1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiLx[1]+2662560.0*rdxCp2R4[0]*phiC[1]-1286983.032055978*rdxCp2R4[0]*bcVals[1]+643491.516027989*phiLx[0]*rdxCp2R4[0])*omega-12960.0*phiC[1]*rdxCp2R4[1]-535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-2662560.0*rdxCp2R4[0]*phiC[1]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiLx[1]+224473.7846609264*rdxCp2[0]*bcVals[1]+((-98207.28078915528*phiUy[0])-14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiLx[1]+2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiLx[0]-2061183.762277152*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiLx[1]+6590799.732961089*rdxCp2R3[0]*bcVals[1]+(3407636.758811008*phiLx[0]-6703036.625291553*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiLx[1]-108000.0*rdxCp2[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiLx[1]-332100.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiLx[1]+241200.0*rdxCp2R3[0]*bcVals[1]+(266400.0*phiLx[0]-387000.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(1052640.0*phiUy[0]-2662560.0*phiC[0])*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiLx[1]+1651200.0*rdxCp2[0]*bcVals[1]+(2563956.0*phiUy[0]+928800.0*phiLx[0]-7679628.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiLx[1]+5004704.0*rdxCp2Sq[0]*bcVals[1]+(532953.0*phiUy[0]+1842246.0*phiLx[0]-3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiLx[1]+1303392.0*rdxCp2R3[0]*bcVals[1]+((-11340.0*phiUy[0])+456408.0*phiLx[0]-535788.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiLx[1]+37440.0*rdxCp2R4[0]*bcVals[1]+(12960.0*phiLx[0]-12960.0*phiC[0])*rdxCp2R4[0])*omega+2662560.0*phiC[0]*rdxCp2R4[1]+7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+12960.0*phiC[0]*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(5263200.0*phiUy[1]-1.33128e+7*phiC[1])*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiLx[1]-3.839814e+7*rdxCp2[0]*phiC[1]+1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(3407636.758811008*phiUy[0]-6703036.625291553*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiLx[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]+1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiUy[0]-2061183.762277152*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiLx[1]-2678940.0*rdxCp2R3[0]*phiC[1]+2731097.713374604*rdxCp2R3[0]*bcVals[1]+((-14029.6115413079*phiUy[0])-98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]+74824.59488697546*rdxCp2R4[0]*bcVals[1])*omega+1.33128e+7*phiC[1]*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiLx[1]+285996.229345773*rdxCp2[0]*bcVals[1]+(1745283.675738703*phiUy[0]+160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiLx[1]+871845.0944978701*rdxCp2Sq[0]*bcVals[1]+(659547.6270141526*phiUy[0]+301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiLx[1]+201610.7140010173*rdxCp2R3[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiLx[1]+688000.0*rdxCp2[0]*bcVals[1]+(266400.0*phiUy[0]-387000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+580800.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiLx[1]+99600.0*rdxCp2R3[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiC[0]-1440.0*phiUy[0])*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiLx[1]-3360.0*rdxCp2[0]*bcVals[1]+((-11799.0*phiUy[0])-1890.0*phiLx[0]+13689.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiLx[1]-33726.0*rdxCp2Sq[0]*bcVals[1]+((-1890.0*phiUy[0])-11799.0*phiLx[0]+13689.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiLx[1]-4160.0*rdxCp2R3[0]*bcVals[1]+(1440.0*phiC[0]-1440.0*phiLx[0])*rdxCp2R3[0])*omega-1440.0*phiC[0]*rdxCp2R3[1]-13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-1440.0*phiC[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiC[1]-1440.0*phiUy[1])*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiLx[1]+13689.0*rdxCp2[0]*phiC[1]-4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiLx[1]+13689.0*rdxCp2Sq[0]*phiC[1]-14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1]-1662.768775266122*rdxCp2R3[0]*bcVals[1])*omega-1440.0*phiC[1]*rdxCp2R3[1]-13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiLx[1]-831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiLx[1]-8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiLx[1]-2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiLx[1]-4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiUy[0]-1125.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(328950.0*phiLy[0]-832050.0*phiC[0])*rdxCp2R3[1]+((-1533436.541464953*rdxCp2[0]*phiLy[1])-335151.8312645776*rdxCp2[0]*phiLx[1]+3096000.0*rdxCp2[0]*bcVals[1]+(2715915.0*phiLy[0]-193500.0*phiLx[0]-8611395.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(90490.99444143593*rdxCp2Sq[0]*phiLy[1]+2176434.423012786*rdxCp2Sq[0]*phiLx[1]+9081960.0*rdxCp2Sq[0]*bcVals[1]+((-193500.0*phiLy[0])+2715915.0*phiLx[0]-8611395.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiLx[1]+1006200.0*rdxCp2R3[0]*bcVals[1]+(328950.0*phiLx[0]-832050.0*phiC[0])*rdxCp2R3[0])*omega+832050.0*phiC[0]*rdxCp2R3[1]+8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+832050.0*phiC[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(832050.0*phiC[1]-328950.0*phiLy[1])*rdxCp2R3[1]+((-125505.0*rdxCp2[0]*phiLy[1])+1290000.0*rdxCp2[0]*phiLx[1]+8611395.0*rdxCp2[0]*phiC[1]-4468691.083527703*rdxCp2[0]*bcVals[1]+(567939.4598018348*phiLy[0]+1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(40635.0*rdxCp2Sq[0]*phiLy[1]+2054550.0*rdxCp2Sq[0]*phiLx[1]+8611395.0*rdxCp2Sq[0]*phiC[1]-4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(1919718.512568964*phiLx[0]-33515.18312645776*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiLx[1]+832050.0*rdxCp2R3[0]*phiC[1]-402182.1975174932*rdxCp2R3[0]*bcVals[1]+201091.0987587466*phiLx[0]*rdxCp2R3[0])*omega-832050.0*phiC[1]*rdxCp2R3[1]-8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-832050.0*rdxCp2R3[0]*phiC[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2]+201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-359640.0*rdxCp2[0]*phiLy[1])-58050.0*rdxCp2[0]*phiLx[1]+536242.9300233242*rdxCp2[0]*bcVals[1]+(1919718.512568964*phiLy[0]-33515.18312645776*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-522450.0*rdxCp2Sq[0]*phiLy[1])+493650.0*rdxCp2Sq[0]*phiLx[1]+1098466.622160182*rdxCp2Sq[0]*bcVals[1]+(1117172.770881926*phiLy[0]+567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+((-437394.7904353685*rdxCp2[0]*phiLy[1])+74478.18472546172*rdxCp2[0]*phiLx[1]-258000.0*rdxCp2[0]*bcVals[1]+(44400.0*phiLy[0]+64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-78202.0939617348*rdxCp2Sq[0]*phiLy[1])+40789.79651824706*rdxCp2Sq[0]*phiLx[1]+40200.0*rdxCp2Sq[0]*bcVals[1]+(64500.0*phiLy[0]+44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(12960.0*phiC[0]-12960.0*phiLy[0])*rdxCp2R4[1]+(176773.1054204795*rdxCp2[0]*phiLy[1]+19641.45615783106*rdxCp2[0]*phiLx[1]-181440.0*rdxCp2[0]*bcVals[1]+((-456408.0*phiLy[0])+11340.0*phiLx[0]+535788.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(814839.8383191626*rdxCp2Sq[0]*phiLy[1]-386469.0325912283*rdxCp2Sq[0]*phiLx[1]-2626452.0*rdxCp2Sq[0]*bcVals[1]+((-1842246.0*phiLy[0])-532953.0*phiLx[0]+3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(434356.7733188925*rdxCp2R3[0]*phiLy[1]-2064285.865273508*rdxCp2R3[0]*phiLx[1]-8373744.0*rdxCp2R3[0]*bcVals[1]+((-928800.0*phiLy[0])-2563956.0*phiLx[0]+7679628.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiLx[1]-3219840.0*rdxCp2R4[0]*bcVals[1]+(2662560.0*phiC[0]-1052640.0*phiLx[0])*rdxCp2R4[0])*omega-12960.0*phiC[0]*rdxCp2R4[1]-535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-2662560.0*phiC[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+((12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(12960.0*phiLy[1]-12960.0*phiC[1])*rdxCp2R4[1]+(157788.0*rdxCp2[0]*phiLy[1]-75600.0*rdxCp2[0]*phiLx[1]-535788.0*rdxCp2[0]*phiC[1]+261886.0821044141*rdxCp2[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(465750.0*rdxCp2Sq[0]*phiLy[1]-727155.0*rdxCp2Sq[0]*phiLx[1]-3688425.0*rdxCp2Sq[0]*phiC[1]+1922680.319449907*rdxCp2Sq[0]*bcVals[1]+((-301792.5327108011*phiLy[0])-659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(195048.0*rdxCp2R3[0]*phiLy[1]-1862820.0*rdxCp2R3[0]*phiLx[1]-7679628.0*rdxCp2R3[0]*phiC[1]+3812313.1094914*rdxCp2R3[0]*bcVals[1]+((-160872.8790069972*phiLy[0])-1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiLx[1]-2662560.0*rdxCp2R4[0]*phiC[1]+1286983.032055978*rdxCp2R4[0]*bcVals[1]-643491.516027989*phiLx[0]*rdxCp2R4[0])*omega+12960.0*phiC[1]*rdxCp2R4[1]+535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+2662560.0*rdxCp2R4[0]*phiC[1])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2]+(218700.0*rdxCp2[0]*phiLy[1]+24300.0*rdxCp2[0]*phiLx[1]-224473.7846609264*rdxCp2[0]*bcVals[1]+(98207.28078915528*phiLy[0]+14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-664200.0*rdxCp2Sq[0]*phiLx[1])-2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLy[0]-814886.6036909672*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3134700.0*rdxCp2R3[0]*phiLy[1])-2961900.0*rdxCp2R3[0]*phiLx[1]-6590799.732961089*rdxCp2R3[0]*bcVals[1]+(6703036.625291553*phiLy[0]-3407636.758811008*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-155884.5726811989*rdxCp2[0]*phiLy[1])-31176.91453623978*rdxCp2[0]*phiLx[1]+108000.0*rdxCp2[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-687061.2540923841*rdxCp2Sq[0]*phiLy[1])-175759.8556980518*rdxCp2Sq[0]*phiLx[1]+332100.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-469212.5637704087*rdxCp2R3[0]*phiLy[1])-244738.7791094823*rdxCp2R3[0]*phiLx[1]-241200.0*rdxCp2R3[0]*bcVals[1]+(387000.0*phiLy[0]-266400.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(2662560.0*phiC[0]-1052640.0*phiLy[0])*rdxCp2R4[1]+((-454351.5678414679*rdxCp2[0]*phiLy[1])-893738.2167055405*rdxCp2[0]*phiLx[1]-1651200.0*rdxCp2[0]*bcVals[1]+((-2563956.0*phiLy[0])-928800.0*phiLx[0]+7679628.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-108651.5471587956*rdxCp2Sq[0]*phiLy[1])-1772702.040022519*rdxCp2Sq[0]*phiLx[1]-5004704.0*rdxCp2Sq[0]*bcVals[1]+((-532953.0*phiLy[0])-1842246.0*phiLx[0]+3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(1870.614872174387*rdxCp2R3[0]*phiLy[1]-439178.8027671645*rdxCp2R3[0]*phiLx[1]-1303392.0*rdxCp2R3[0]*bcVals[1]+(11340.0*phiLy[0]-456408.0*phiLx[0]+535788.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiLx[1]-37440.0*rdxCp2R4[0]*bcVals[1]+(12960.0*phiC[0]-12960.0*phiLx[0])*rdxCp2R4[0])*omega-2662560.0*phiC[0]*rdxCp2R4[1]-7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-12960.0*phiC[0]*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(1.33128e+7*phiC[1]-5263200.0*phiLy[1])*rdxCp2R4[1]+((-1408860.0*rdxCp2[0]*phiLy[1])+6450000.0*rdxCp2[0]*phiLx[1]+3.839814e+7*rdxCp2[0]*phiC[1]-1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(6703036.625291553*phiLx[0]-3407636.758811008*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+(63990.0*rdxCp2Sq[0]*phiLy[1]+1983375.0*rdxCp2Sq[0]*phiLx[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]-1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLx[0]-814886.6036909672*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(9720.0*rdxCp2R3[0]*phiLy[1]+94500.0*rdxCp2R3[0]*phiLx[1]+2678940.0*rdxCp2R3[0]*phiC[1]-2731097.713374604*rdxCp2R3[0]*bcVals[1]+(14029.6115413079*phiLy[0]+98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1]-74824.59488697546*rdxCp2R4[0]*bcVals[1])*omega-1.33128e+7*phiC[1]*rdxCp2R4[1]-3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]-643491.516027989*phiLy[0]*rdxCp2R4[1]+((-106560.0*rdxCp2[0]*phiLy[1])-154800.0*rdxCp2[0]*phiLx[1]-285996.229345773*rdxCp2[0]*bcVals[1]+((-1745283.675738703*phiLy[0])-160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-66420.0*rdxCp2Sq[0]*phiLy[1])-290400.0*rdxCp2Sq[0]*phiLx[1]-871845.0944978701*rdxCp2Sq[0]*bcVals[1]+((-659547.6270141526*phiLy[0])-301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-10800.0*rdxCp2R3[0]*phiLy[1])-63000.0*rdxCp2R3[0]*phiLx[1]-201610.7140010173*rdxCp2R3[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+((-2016730.678300898*rdxCp2[0]*phiLy[1])+372390.9236273086*rdxCp2[0]*phiLx[1]-688000.0*rdxCp2[0]*bcVals[1]+(387000.0*phiLx[0]-266400.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-543205.7742697513*rdxCp2Sq[0]*phiLy[1])-580800.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-18706.14872174387*rdxCp2R3[0]*phiLy[1])-25980.76211353316*rdxCp2R3[0]*phiLx[1]-99600.0*rdxCp2R3[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(1440.0*phiLy[0]-1440.0*phiC[0])*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiLy[1]+1818.653347947321*rdxCp2[0]*phiLx[1]+3360.0*rdxCp2[0]*bcVals[1]+(11799.0*phiLy[0]+1890.0*phiLx[0]-13689.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2Sq[0]*phiLy[1]+11353.59304361399*rdxCp2Sq[0]*phiLx[1]+33726.0*rdxCp2Sq[0]*bcVals[1]+(1890.0*phiLy[0]+11799.0*phiLx[0]-13689.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiLx[1]+4160.0*rdxCp2R3[0]*bcVals[1]+(1440.0*phiLx[0]-1440.0*phiC[0])*rdxCp2R3[0])*omega+1440.0*phiC[0]*rdxCp2R3[1]+13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+1440.0*phiC[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(1440.0*phiLy[1]-1440.0*phiC[1])*rdxCp2R3[1]+(2664.0*rdxCp2[0]*phiLy[1]-2625.0*rdxCp2[0]*phiLx[1]-13689.0*rdxCp2[0]*phiC[1]+4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLy[0]-2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(324.0*rdxCp2Sq[0]*phiLy[1]-450.0*rdxCp2Sq[0]*phiLx[1]-13689.0*rdxCp2Sq[0]*phiC[1]+14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLy[0]-467.6537180435967*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]+1662.768775266122*rdxCp2R3[0]*bcVals[1])*omega+1440.0*phiC[1]*rdxCp2R3[1]+13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]+(450.0*rdxCp2[0]*phiLy[1]+450.0*rdxCp2[0]*phiLx[1]+831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiLy[1])+2625.0*rdxCp2Sq[0]*phiLx[1]+8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-4546.633369868302*rdxCp2[0]*phiLy[1])-1082.531754730548*rdxCp2[0]*phiLx[1]+2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLy[0]-1125.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-779.4228634059946*rdxCp2Sq[0]*phiLy[1])+1082.531754730548*rdxCp2Sq[0]*phiLx[1]+4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 
  double *dxUy = dx[3]; 
  double *dxLy = dx[4]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdxCp2[2]; 
  double rdxLx[2]; 
  double rdxUx[2]; 
  double rdxLxSq[2]; 
  double rdxUxSq[2]; 
  double rdxLxCu[2]; 
  double rdxUxCu[2]; 
  double rdxLxR4[2]; 
  double rdxUxR4[2]; 
  double rdxLy[2]; 
  double rdxUy[2]; 
  double rdxLySq[2]; 
  double rdxUySq[2]; 
  double rdxLyCu[2]; 
  double rdxUyCu[2]; 
  double rdxLyR4[2]; 
  double rdxUyR4[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = volFac*4.0/(dxLx[1]*dxLx[1]); 
  rdxUx[1]   = volFac*4.0/(dxUx[1]*dxUx[1]); 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = volFac*4.0/(dxLy[0]*dxLy[0]); 
  rdxUy[0]   = volFac*4.0/(dxUy[0]*dxUy[0]); 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = volFac*4.0/(dxLy[1]*dxLy[1]); 
  rdxUy[1]   = volFac*4.0/(dxUy[1]*dxUy[1]); 
  rdxLySq[1] = rdxLy[1]*rdxLy[1]; 
  rdxUySq[1] = rdxUy[1]*rdxUy[1]; 
  rdxLyCu[1] = rdxLySq[1]*rdxLy[1]; 
  rdxUyCu[1] = rdxUySq[1]*rdxUy[1]; 
  rdxLyR4[1] = rdxLyCu[1]*rdxLy[1]; 
  rdxUyR4[1] = rdxUyCu[1]*rdxUy[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((((9600.0*rdxUx[0]-9600.0*rdxLx[0])*rdxUySq[1]+(9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUy[1]+(9600.0*rdxLx[0]-9600.0*rdxUx[0])*rdxLySq[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLy[1])*rho[3]+((-415.6921938165305*rdxUyCu[1])+((-27435.68479189101*rdxLy[1])-25495.78788741387*rdxUx[0]-25495.78788741387*rdxLx[0])*rdxUySq[1]+(27435.68479189101*rdxLySq[1]-25080.09569359734*rdxUxSq[0]-23140.1987891202*rdxLx[0]*rdxUx[0]-25080.09569359734*rdxLxSq[0])*rdxUy[1]+415.6921938165305*rdxLyCu[1]+(25495.78788741387*rdxUx[0]+25495.78788741387*rdxLx[0])*rdxLySq[1]+(25080.09569359734*rdxUxSq[0]+23140.1987891202*rdxLx[0]*rdxUx[0]+25080.09569359734*rdxLxSq[0])*rdxLy[1])*rho[2]+((25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxUySq[1]+((23140.1987891202*rdxLx[0]-23140.1987891202*rdxUx[0])*rdxLy[1]-25495.78788741387*rdxUxSq[0]+25495.78788741387*rdxLxSq[0])*rdxUy[1]+(25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxLySq[1]+(25495.78788741387*rdxLxSq[0]-25495.78788741387*rdxUxSq[0])*rdxLy[1]-415.6921938165305*rdxUxCu[0]-27435.68479189101*rdxLx[0]*rdxUxSq[0]+27435.68479189101*rdxLxSq[0]*rdxUx[0]+415.6921938165305*rdxLxCu[0])*rho[1]+1104.0*rho[0]*rdxUyCu[1]+(75072.0*rho[0]*rdxLy[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0])*rdxUySq[1]+(75072.0*rho[0]*rdxLySq[1]+(164368.0*rdxUx[0]+164368.0*rdxLx[0])*rho[0]*rdxLy[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0])*rdxUy[1]+1104.0*rho[0]*rdxLyCu[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0]*rdxLySq[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0]*rdxLy[1]+(1104.0*rdxUxCu[0]+75072.0*rdxLx[0]*rdxUxSq[0]+75072.0*rdxLxSq[0]*rdxUx[0]+1104.0*rdxLxCu[0])*rho[0])*volFac+((9375.0*rdxUx[0]-9375.0*rdxLx[0])*rdxUyCu[1]+((12525.0*rdxUx[0]-12525.0*rdxLx[0])*rdxLy[1]+9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+((17775.0*rdxUx[0]-17775.0*rdxLx[0])*rdxLySq[1]+(18000.0*rdxUxSq[0]-18000.0*rdxLxSq[0])*rdxLy[1]+225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]+9600.0*rdxUxSq[0]+18000.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])+9375.0*rdxUxCu[0]+12525.0*rdxLx[0]*rdxUxSq[0]+17775.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+((-9600.0*rdxUxSq[0])-18000.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-9375.0*rdxUxCu[0])-12525.0*rdxLx[0]*rdxUxSq[0]-17775.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((17775.0*rdxLx[0]-17775.0*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((12525.0*rdxLx[0]-12525.0*rdxUx[0])*rdxLySq[1]+(18000.0*rdxLxSq[0]-18000.0*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9375.0*rdxLx[0]-9375.0*rdxUx[0])*rdxLyCu[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLySq[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-18000.0*rdxLx[0]*rdxUx[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-17775.0*rdxLx[0]*rdxUxSq[0]-12525.0*rdxLxSq[0]*rdxUx[0]-9375.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(18000.0*rdxLx[0]*rdxUx[0]+9600.0*rdxLxSq[0])*rdxLySq[1]+(17775.0*rdxLx[0]*rdxUxSq[0]+12525.0*rdxLxSq[0]*rdxUx[0]+9375.0*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((-415.6921938165305*rdxUyR4[1])+((-28630.79984911354*rdxLy[1])-25729.61474643567*rdxUx[0]-25729.61474643567*rdxLx[0])*rdxUyCu[1]+((-52637.02404201817*rdxLySq[1])+((-88966.78973077537*rdxUx[0])-88966.78973077537*rdxLx[0])*rdxLy[1]-25911.4800812304*rdxUxSq[0]-78842.95276053529*rdxLx[0]*rdxUx[0]-25911.4800812304*rdxLxSq[0])*rdxUySq[1]+((-779.4228634059946*rdxLyCu[1])+((-48038.4291479228*rdxUx[0])-48038.4291479228*rdxLx[0])*rdxLySq[1]+((-47856.56381312806*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0]-47856.56381312806*rdxLxSq[0])*rdxLy[1]-597.5575286112626*rdxUxCu[0]-40633.91194556586*rdxLx[0]*rdxUxSq[0]-40633.91194556586*rdxLxSq[0]*rdxUx[0]-597.5575286112626*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-233.8268590217983*rdxUx[0]*rdxUyCu[1])+((-15432.57269543869*rdxUx[0]*rdxLy[1])-9145.22826396367*rdxUxSq[0]-19537.53310937693*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15432.57269543869*rdxUx[0]*rdxLySq[1]-8911.401404941873*rdxUxCu[0]-13016.36181888011*rdxLx[0]*rdxUxSq[0]-19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]+19537.53310937693*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8911.401404941873*rdxUxCu[0]+13016.36181888011*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+(779.4228634059946*rdxLy[1]*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(48038.4291479228*rdxUx[0]+48038.4291479228*rdxLx[0])*rdxLy[1])*rdxUySq[1]+(28630.79984911354*rdxLyCu[1]+(88966.78973077537*rdxUx[0]+88966.78973077537*rdxLx[0])*rdxLySq[1]+(47856.56381312806*rdxUxSq[0]+99090.62670101546*rdxLx[0]*rdxUx[0]+47856.56381312806*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+415.6921938165305*rdxLyR4[1]+(25729.61474643567*rdxUx[0]+25729.61474643567*rdxLx[0])*rdxLyCu[1]+(25911.4800812304*rdxUxSq[0]+78842.95276053529*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*rdxLySq[1]+(597.5575286112626*rdxUxCu[0]+40633.91194556586*rdxLx[0]*rdxUxSq[0]+40633.91194556586*rdxLxSq[0]*rdxUx[0]+597.5575286112626*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15432.57269543869*rdxLx[0]*rdxLy[1])-19537.53310937693*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(15432.57269543869*rdxLx[0]*rdxLySq[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-13016.36181888011*rdxLxSq[0]*rdxUx[0]-8911.401404941873*rdxLxCu[0])*rdxUy[1]+233.8268590217983*rdxLx[0]*rdxLyCu[1]+(19537.53310937693*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(19303.70625035514*rdxLx[0]*rdxUxSq[0]+13016.36181888011*rdxLxSq[0]*rdxUx[0]+8911.401404941873*rdxLxCu[0])*rdxLy[1])*phiLx[2]+396.0*phiUy[0]*rdxUyR4[1]+((27378.0*phiUy[0]+846.0*phiLy[0])*rdxLy[1]+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiUy[1]-597.5575286112626*rdxUx[0]*phiUx[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(24531.0*phiUy[0]+621.0*phiUx[0])*rdxUx[0]+(24531.0*phiUy[0]+621.0*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((57078.0*phiUy[0]+57078.0*phiLy[0])*rdxLySq[1]+((13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(92457.0*phiUy[0]+42228.0*phiUx[0]+52131.0*phiLy[0])*rdxUx[0]+(92457.0*phiUy[0]+52131.0*phiLy[0]+42228.0*phiLx[0])*rdxLx[0])*rdxLy[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiUy[1]+((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUy[0]+24756.0*phiUx[0])*rdxUxSq[0]+(79932.0*phiUy[0]+51906.0*phiUx[0]+51906.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiUy[0]+24756.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((846.0*phiUy[0]+27378.0*phiLy[0])*rdxLyCu[1]+((19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(52131.0*phiUy[0]+42228.0*phiUx[0]+92457.0*phiLy[0])*rdxUx[0]+(52131.0*phiUy[0]+92457.0*phiLy[0]+42228.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiUy[1]+((-78842.95276053529*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0])*phiUx[1]+(19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiLy[1]+(99090.62670101546*rdxLx[0]*rdxUx[0]+78842.95276053529*rdxLxSq[0])*phiLx[1]+(51906.0*phiUy[0]+79932.0*phiUx[0]+51906.0*phiLy[0])*rdxUxSq[0]+(104982.0*phiUy[0]+104982.0*phiUx[0]+104982.0*phiLy[0]+104982.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(51906.0*phiUy[0]+51906.0*phiLy[0]+79932.0*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiUy[1]+((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(621.0*phiUy[0]+24531.0*phiUx[0])*rdxUxCu[0]+(42228.0*phiUy[0]+92457.0*phiUx[0]+52131.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(42228.0*phiUy[0]+52131.0*phiUx[0]+92457.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiUy[0]+24531.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+396.0*phiLy[0]*rdxLyR4[1]+((-597.5575286112626*rdxUx[0]*phiUx[1])+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiLy[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(621.0*phiUx[0]+24531.0*phiLy[0])*rdxUx[0]+(24531.0*phiLy[0]+621.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+(((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiLy[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUx[0]+24756.0*phiLy[0])*rdxUxSq[0]+(51906.0*phiUx[0]+79932.0*phiLy[0]+51906.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiLy[0]+24756.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiLy[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(24531.0*phiUx[0]+621.0*phiLy[0])*rdxUxCu[0]+(92457.0*phiUx[0]+42228.0*phiLy[0]+52131.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(52131.0*phiUx[0]+42228.0*phiLy[0]+92457.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiLy[0]+24531.0*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-415.6921938165305*rdxUxR4[0])-28630.79984911354*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(779.4228634059946*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+28630.79984911354*rdxLxCu[0]*rdxUx[0]+415.6921938165305*rdxLxR4[0])*phiLx[1]+396.0*phiUx[0]*rdxUxR4[0]+(27378.0*phiUx[0]+846.0*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(57078.0*phiUx[0]+57078.0*phiLx[0])*rdxLxSq[0]*rdxUxSq[0]+(846.0*phiUx[0]+27378.0*phiLx[0])*rdxLxCu[0]*rdxUx[0]+396.0*phiLx[0]*rdxLxR4[0])/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxUyCu[1]+(27435.68479189101*rdxLy[1]+9976.61265159673*rdxUx[0]+9976.61265159673*rdxLx[0])*rdxUySq[1]+((-27435.68479189101*rdxLySq[1])+9560.920457780201*rdxUxSq[0]-7898.15168251408*rdxLx[0]*rdxUx[0]+9560.920457780201*rdxLxSq[0])*rdxUy[1]-415.6921938165305*rdxLyCu[1]+((-9976.61265159673*rdxUx[0])-9976.61265159673*rdxLx[0])*rdxLySq[1]+((-9560.920457780201*rdxUxSq[0])+7898.15168251408*rdxLx[0]*rdxUx[0]-9560.920457780201*rdxLxSq[0])*rdxLy[1])*rho[3]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[2]+((-1104.0*rdxUyCu[1])+((-75072.0*rdxLy[1])-27600.0*rdxUx[0]-27600.0*rdxLx[0])*rdxUySq[1]+((-75072.0*rdxLySq[1])+((-126960.0*rdxUx[0])-126960.0*rdxLx[0])*rdxLy[1]-26928.0*rdxUxSq[0]-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxUy[1]-1104.0*rdxLyCu[1]+((-27600.0*rdxUx[0])-27600.0*rdxLx[0])*rdxLySq[1]+((-26928.0*rdxUxSq[0])-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxLy[1]-432.0*rdxUxCu[0]-29376.0*rdxLx[0]*rdxUxSq[0]-29376.0*rdxLxSq[0]*rdxUx[0]-432.0*rdxLxCu[0])*rho[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxUySq[1]+((60164.51685171252*rdxUx[0]-60164.51685171252*rdxLx[0])*rho[0]*rdxLy[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0])*rdxUy[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxLySq[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0]*rdxLy[1]+(1080.799703922979*rdxUxCu[0]+71332.78045891662*rdxLx[0]*rdxUxSq[0]-71332.78045891662*rdxLxSq[0]*rdxUx[0]-1080.799703922979*rdxLxCu[0])*rho[0])*volFac+(415.6921938165305*rdxUyR4[1]+(28630.79984911354*rdxLy[1]+10574.170180208*rdxUx[0]+10574.170180208*rdxLx[0])*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(68719.1157902952*rdxUx[0]+68719.1157902952*rdxLx[0])*rdxLy[1]+10392.30484541326*rdxUxSq[0]+47804.60228890101*rdxLx[0]*rdxUx[0]+10392.30484541326*rdxLxSq[0])*rdxUySq[1]+(779.4228634059946*rdxLyCu[1]+(19303.70625035514*rdxUx[0]+19303.70625035514*rdxLx[0])*rdxLySq[1]+(18758.11024597094*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0]+18758.11024597094*rdxLxSq[0])*rdxLy[1]+233.8268590217983*rdxUxCu[0]+15900.22641348229*rdxLx[0]*rdxUxSq[0]+15900.22641348229*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-181.8653347947321*rdxUx[0]*rdxUyCu[1])+((-12003.11209645232*rdxUx[0]*rdxLy[1])+9145.22826396367*rdxUxSq[0]-17874.76433411081*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(12003.11209645232*rdxUx[0]*rdxLySq[1]+9327.093598758403*rdxUxCu[0]+3455.441361099909*rdxLx[0]*rdxUxSq[0]-17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+181.8653347947321*rdxUx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxUxSq[0])*rdxLySq[1]+((-9327.093598758403*rdxUxCu[0])-3455.441361099909*rdxLx[0]*rdxUxSq[0]+17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((-779.4228634059946*rdxLy[1]*rdxUyCu[1])+(((-19303.70625035514*rdxUx[0])-19303.70625035514*rdxLx[0])*rdxLy[1]-52637.02404201817*rdxLySq[1])*rdxUySq[1]+((-28630.79984911354*rdxLyCu[1])+((-68719.1157902952*rdxUx[0])-68719.1157902952*rdxLx[0])*rdxLySq[1]+((-18758.11024597094*rdxUxSq[0])-40893.71956670118*rdxLx[0]*rdxUx[0]-18758.11024597094*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-415.6921938165305*rdxLyR4[1]+((-10574.170180208*rdxUx[0])-10574.170180208*rdxLx[0])*rdxLyCu[1]+((-10392.30484541326*rdxUxSq[0])-47804.60228890101*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-233.8268590217983*rdxUxCu[0])-15900.22641348229*rdxLx[0]*rdxUxSq[0]-15900.22641348229*rdxLxSq[0]*rdxUx[0]-233.8268590217983*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-181.8653347947321*rdxLx[0]*rdxUyCu[1])+((-12003.11209645232*rdxLx[0]*rdxLy[1])-17874.76433411081*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(12003.11209645232*rdxLx[0]*rdxLySq[1]-17692.89899931608*rdxLx[0]*rdxUxSq[0]+3455.441361099909*rdxLxSq[0]*rdxUx[0]+9327.093598758403*rdxLxCu[0])*rdxUy[1]+181.8653347947321*rdxLx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(17692.89899931608*rdxLx[0]*rdxUxSq[0]-3455.441361099909*rdxLxSq[0]*rdxUx[0]-9327.093598758403*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((24375.0*rdxLx[0]-24375.0*rdxUx[0])*rdxUyCu[1]+((32565.0*rdxLx[0]-32565.0*rdxUx[0])*rdxLy[1]-24960.0*rdxUxSq[0]+24960.0*rdxLxSq[0])*rdxUySq[1]+((46215.0*rdxLx[0]-46215.0*rdxUx[0])*rdxLySq[1]+(46800.0*rdxLxSq[0]-46800.0*rdxUxSq[0])*rdxLy[1]-585.0*rdxUxCu[0]-38610.0*rdxLx[0]*rdxUxSq[0]+38610.0*rdxLxSq[0]*rdxUx[0]+585.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]-8640.0*rdxUxSq[0]+19440.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])-8865.0*rdxUxCu[0]-4275.0*rdxLx[0]*rdxUxSq[0]+19215.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+(8640.0*rdxUxSq[0]-19440.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8865.0*rdxUxCu[0]+4275.0*rdxLx[0]*rdxUxSq[0]-19215.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+((46215.0*rdxUx[0]-46215.0*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((32565.0*rdxUx[0]-32565.0*rdxLx[0])*rdxLySq[1]+(46800.0*rdxUxSq[0]-46800.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(24375.0*rdxUx[0]-24375.0*rdxLx[0])*rdxLyCu[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLySq[1]+(585.0*rdxUxCu[0]+38610.0*rdxLx[0]*rdxUxSq[0]-38610.0*rdxLxSq[0]*rdxUx[0]-585.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-19440.0*rdxLx[0]*rdxUx[0]+8640.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-19215.0*rdxLx[0]*rdxUxSq[0]+4275.0*rdxLxSq[0]*rdxUx[0]+8865.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(19440.0*rdxLx[0]*rdxUx[0]-8640.0*rdxLxSq[0])*rdxLySq[1]+(19215.0*rdxLx[0]*rdxUxSq[0]-4275.0*rdxLxSq[0]*rdxUx[0]-8865.0*rdxLxCu[0])*rdxLy[1])*phiLx[2]-396.0*phiUy[1]*rdxUyR4[1]+(((-27378.0*phiUy[1])-846.0*phiLy[1])*rdxLy[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiUy[1]+483.0*rdxUx[0]*phiUx[1]+483.0*rdxLx[0]*phiLx[1]+(23169.64365284887*phiUy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiUy[0])*rdxLx[0])*rdxUyCu[1]+(((-57078.0*phiUy[1])-57078.0*phiLy[1])*rdxLySq[1]+(((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(33842.54072908828*phiUy[0]-40633.91194556586*phiUx[0]+50189.63625092335*phiLy[0])*rdxUx[0]+((-33842.54072908828*phiUy[0])-50189.63625092335*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLy[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiUy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(23777.59348630554*phiUy[0]+21740.70173660454*phiUx[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiUy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-846.0*phiUy[1])-27378.0*phiLy[1])*rdxLyCu[1]+(((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(50189.63625092335*phiUy[0]-40633.91194556586*phiUx[0]+33842.54072908828*phiLy[0])*rdxUx[0]+((-50189.63625092335*phiUy[0])-33842.54072908828*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLySq[1]+(((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiUy[1]+(22980.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0])*phiUx[1]+((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiLy[1]+(88110.0*rdxLx[0]*rdxUx[0]+22980.0*rdxLxSq[0])*phiLx[1]+(50797.58608438003*phiUy[0]-34876.57506120691*phiUx[0]+50797.58608438003*phiLy[0])*rdxUxSq[0]+(102561.6565193835*phiLx[0]-102561.6565193835*phiUx[0])*rdxLx[0]*rdxUx[0]+((-50797.58608438003*phiUy[0])-50797.58608438003*phiLy[0]+34876.57506120691*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiUy[1]+((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(607.9498334566757*phiUy[0]+22712.38223965068*phiUx[0])*rdxUxCu[0]+(40124.68900814059*phiUy[0]-44349.16092780109*phiUx[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40124.68900814059*phiUy[0])-51862.79733103487*phiUx[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiUy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxUy[1]-396.0*phiLy[1]*rdxLyR4[1]+(483.0*rdxUx[0]*phiUx[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiLy[1]+483.0*rdxLx[0]*phiLx[1]+(23169.64365284887*phiLy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiLy[0])*rdxLx[0])*rdxLyCu[1]+((47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiLy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(21740.70173660454*phiUx[0]+23777.59348630554*phiLy[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiLy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiLy[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(22712.38223965068*phiUx[0]+607.9498334566757*phiLy[0])*rdxUxCu[0]+((-44349.16092780109*phiUx[0])+40124.68900814059*phiLy[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-51862.79733103487*phiUx[0])-40124.68900814059*phiLy[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiLy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-396.0*rdxUxR4[0])-25758.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]+774.0*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(774.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]-25758.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[1]+374.1229744348773*phiUx[0]*rdxUxR4[0]+(24224.46259465831*phiUx[0]+841.7766924784738*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(56024.91542162288*phiLx[0]-56024.91542162288*phiUx[0])*rdxLxSq[0]*rdxUxSq[0]+((-841.7766924784738*phiUx[0])-24224.46259465831*phiLx[0])*rdxLxCu[0]*rdxUx[0]-374.1229744348773*phiLx[0]*rdxLxR4[0]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[2] = -(1.0*((((9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxUySq[1]+((7898.15168251408*rdxLx[0]-7898.15168251408*rdxUx[0])*rdxLy[1]+9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxUy[1]+(9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxLySq[1]+(9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxLy[1]+415.6921938165305*rdxUxCu[0]+27435.68479189101*rdxLx[0]*rdxUxSq[0]-27435.68479189101*rdxLxSq[0]*rdxUx[0]-415.6921938165305*rdxLxCu[0])*rho[3]+((-432.0*rdxUyCu[1])+((-29376.0*rdxLy[1])-26928.0*rdxUx[0]-26928.0*rdxLx[0])*rdxUySq[1]+((-29376.0*rdxLySq[1])+((-81936.0*rdxUx[0])-81936.0*rdxLx[0])*rdxLy[1]-27600.0*rdxUxSq[0]-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxUy[1]-432.0*rdxLyCu[1]+((-26928.0*rdxUx[0])-26928.0*rdxLx[0])*rdxLySq[1]+((-27600.0*rdxUxSq[0])-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxLy[1]-1104.0*rdxUxCu[0]-75072.0*rdxLx[0]*rdxUxSq[0]-75072.0*rdxLxSq[0]*rdxUx[0]-1104.0*rdxLxCu[0])*rho[2]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[1]+1080.799703922979*rho[0]*rdxUyCu[1]+(71332.78045891662*rho[0]*rdxLy[1]+(66289.04850727606*rdxUx[0]+66289.04850727606*rdxLx[0])*rho[0])*rdxUySq[1]+((65208.24880335309*rdxUxSq[0]+60164.51685171252*rdxLx[0]*rdxUx[0]+65208.24880335309*rdxLxSq[0])*rho[0]-71332.78045891662*rho[0]*rdxLySq[1])*rdxUy[1]-1080.799703922979*rho[0]*rdxLyCu[1]+((-66289.04850727606*rdxUx[0])-66289.04850727606*rdxLx[0])*rho[0]*rdxLySq[1]+((-65208.24880335309*rdxUxSq[0])-60164.51685171252*rdxLx[0]*rdxUx[0]-65208.24880335309*rdxLxSq[0])*rho[0]*rdxLy[1])*volFac+((9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxUyCu[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLy[1]+9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1]-181.8653347947321*rdxUxCu[0]-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(233.8268590217983*rdxUx[0]*rdxUyCu[1]+(15900.22641348229*rdxUx[0]*rdxLy[1]+10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15900.22641348229*rdxUx[0]*rdxLySq[1]+(47804.60228890101*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0])*rdxLy[1]+10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+415.6921938165305*rdxUxR4[0]+28630.79984911354*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+((-181.8653347947321*rdxUxCu[0])-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15900.22641348229*rdxLx[0]*rdxLy[1])-18758.11024597094*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxUySq[1]+((-15900.22641348229*rdxLx[0]*rdxLySq[1])+((-40893.71956670118*rdxLx[0]*rdxUx[0])-47804.60228890101*rdxLxSq[0])*rdxLy[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxUy[1]-233.8268590217983*rdxLx[0]*rdxLyCu[1]+((-18758.11024597094*rdxLx[0]*rdxUx[0])-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-19303.70625035514*rdxLx[0]*rdxUxSq[0])-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxLy[1]-779.4228634059946*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-28630.79984911354*rdxLxCu[0]*rdxUx[0]-415.6921938165305*rdxLxR4[0])*phiLx[3]+((-396.0*rdxUyR4[1])+((-25758.0*rdxLy[1])-24099.0*rdxUx[0]-24099.0*rdxLx[0])*rdxUyCu[1]+(51462.0*rdxLySq[1]+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLy[1]-23220.0*rdxUxSq[0]+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxUySq[1]+(774.0*rdxLyCu[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1]+483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-243.0*rdxUx[0]*rdxUyCu[1])+((-16524.0*rdxUx[0]*rdxLy[1])-9972.0*rdxUxSq[0]-20322.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-16524.0*rdxUx[0]*rdxLySq[1])+((-50364.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0])*rdxLy[1]-10125.0*rdxUxCu[0]-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-243.0*rdxUx[0]*rdxLyCu[1]+((-9972.0*rdxUxSq[0])-20322.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-10125.0*rdxUxCu[0])-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-396.0*rdxUxR4[0]-27378.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-846.0*rdxLxCu[0]*rdxUx[0])*phiUx[2]+(774.0*rdxLy[1]*rdxUyCu[1]+(51462.0*rdxLySq[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLy[1])*rdxUySq[1]+((-25758.0*rdxLyCu[1])+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-396.0*rdxLyR4[1]+((-24099.0*rdxUx[0])-24099.0*rdxLx[0])*rdxLyCu[1]+((-23220.0*rdxUxSq[0])+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxLySq[1]+(483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-243.0*rdxLx[0]*rdxUyCu[1])+((-16524.0*rdxLx[0]*rdxLy[1])-20322.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*rdxUySq[1]+((-16524.0*rdxLx[0]*rdxLySq[1])+((-41814.0*rdxLx[0]*rdxUx[0])-50364.0*rdxLxSq[0])*rdxLy[1]-20925.0*rdxLx[0]*rdxUxSq[0]-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxUy[1]-243.0*rdxLx[0]*rdxLyCu[1]+((-20322.0*rdxLx[0]*rdxUx[0])-9972.0*rdxLxSq[0])*rdxLySq[1]+((-20925.0*rdxLx[0]*rdxUxSq[0])-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxLy[1]-846.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-27378.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[2]+374.1229744348773*phiUy[0]*rdxUyR4[1]+((24224.46259465831*phiUy[0]+841.7766924784738*phiLy[0])*rdxLy[1]+(8865.0*rdxLx[0]-8865.0*rdxUx[0])*phiUy[1]-585.0*rdxUx[0]*phiUx[1]+585.0*rdxLx[0]*phiLx[1]+(22712.38223965068*phiUy[0]+607.9498334566757*phiUx[0])*rdxUx[0]+(22712.38223965068*phiUy[0]+607.9498334566757*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((56024.91542162288*phiLy[0]-56024.91542162288*phiUy[0])*rdxLySq[1]+((4275.0*rdxLx[0]-4275.0*rdxUx[0])*phiUy[1]-38610.0*rdxUx[0]*phiUx[1]+(19215.0*rdxLx[0]-19215.0*rdxUx[0])*phiLy[1]+38610.0*rdxLx[0]*phiLx[1]+((-44349.16092780109*phiUy[0])+40124.68900814059*phiUx[0]+51862.79733103487*phiLy[0])*rdxUx[0]+((-44349.16092780109*phiUy[0])+51862.79733103487*phiLy[0]+40124.68900814059*phiLx[0])*rdxLx[0])*rdxLy[1]+(8640.0*rdxLxSq[0]-8640.0*rdxUxSq[0])*phiUy[1]+((-24960.0*rdxUxSq[0])-46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(46800.0*rdxLx[0]*rdxUx[0]+24960.0*rdxLxSq[0])*phiLx[1]+(21740.70173660454*phiUy[0]+23777.59348630554*phiUx[0])*rdxUxSq[0]+((-34876.57506120691*phiUy[0])+50797.58608438003*phiUx[0]+50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+(21740.70173660454*phiUy[0]+23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-841.7766924784738*phiUy[0])-24224.46259465831*phiLy[0])*rdxLyCu[1]+((19215.0*rdxUx[0]-19215.0*rdxLx[0])*phiUy[1]+38610.0*rdxUx[0]*phiUx[1]+(4275.0*rdxUx[0]-4275.0*rdxLx[0])*phiLy[1]-38610.0*rdxLx[0]*phiLx[1]+((-51862.79733103487*phiUy[0])-40124.68900814059*phiUx[0]+44349.16092780109*phiLy[0])*rdxUx[0]+((-51862.79733103487*phiUy[0])+44349.16092780109*phiLy[0]-40124.68900814059*phiLx[0])*rdxLx[0])*rdxLySq[1]+((19440.0*rdxUxSq[0]-19440.0*rdxLxSq[0])*phiUy[1]+(19440.0*rdxLxSq[0]-19440.0*rdxUxSq[0])*phiLy[1]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxUxSq[0]+(102561.6565193835*phiLy[0]-102561.6565193835*phiUy[0])*rdxLx[0]*rdxUx[0]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxLxSq[0])*rdxLy[1]+(225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*phiUy[1]+((-24375.0*rdxUxCu[0])-32565.0*rdxLx[0]*rdxUxSq[0]-46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(46215.0*rdxLx[0]*rdxUxSq[0]+32565.0*rdxLxSq[0]*rdxUx[0]+24375.0*rdxLxCu[0])*phiLx[1]+(23169.64365284887*phiUx[0]-597.5575286112626*phiUy[0])*rdxUxCu[0]+((-40633.91194556586*phiUy[0])+33842.54072908828*phiUx[0]+50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40633.91194556586*phiUy[0])+50189.63625092335*phiUx[0]+33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(23169.64365284887*phiLx[0]-597.5575286112626*phiUy[0])*rdxLxCu[0])*rdxUy[1]-374.1229744348773*phiLy[0]*rdxLyR4[1]+(585.0*rdxUx[0]*phiUx[1]+(8865.0*rdxUx[0]-8865.0*rdxLx[0])*phiLy[1]-585.0*rdxLx[0]*phiLx[1]+((-607.9498334566757*phiUx[0])-22712.38223965068*phiLy[0])*rdxUx[0]+((-22712.38223965068*phiLy[0])-607.9498334566757*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((24960.0*rdxUxSq[0]+46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(8640.0*rdxUxSq[0]-8640.0*rdxLxSq[0])*phiLy[1]+((-46800.0*rdxLx[0]*rdxUx[0])-24960.0*rdxLxSq[0])*phiLx[1]+((-23777.59348630554*phiUx[0])-21740.70173660454*phiLy[0])*rdxUxSq[0]+((-50797.58608438003*phiUx[0])+34876.57506120691*phiLy[0]-50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+((-21740.70173660454*phiLy[0])-23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((24375.0*rdxUxCu[0]+32565.0*rdxLx[0]*rdxUxSq[0]+46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*phiLy[1]+((-46215.0*rdxLx[0]*rdxUxSq[0])-32565.0*rdxLxSq[0]*rdxUx[0]-24375.0*rdxLxCu[0])*phiLx[1]+(597.5575286112626*phiLy[0]-23169.64365284887*phiUx[0])*rdxUxCu[0]+((-33842.54072908828*phiUx[0])+40633.91194556586*phiLy[0]-50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-50189.63625092335*phiUx[0])+40633.91194556586*phiLy[0]-33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(597.5575286112626*phiLy[0]-23169.64365284887*phiLx[0])*rdxLxCu[0])*rdxLy[1]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[3] = (((144.0*rdxUyCu[1]+(9792.0*rdxLy[1]+3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxUySq[1]+(9792.0*rdxLySq[1]+(31568.0*rdxUx[0]+31568.0*rdxLx[0])*rdxLy[1]+3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxUy[1]+144.0*rdxLyCu[1]+(3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxLySq[1]+(3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxLy[1]+144.0*rdxUxCu[0]+9792.0*rdxLx[0]*rdxUxSq[0]+9792.0*rdxLxSq[0]*rdxUx[0]+144.0*rdxLxCu[0])*rho[3]+((8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxUySq[1]+((6845.064791512203*rdxUx[0]-6845.064791512203*rdxLx[0])*rdxLy[1]-8646.397631383834*rdxUxSq[0]+8646.397631383834*rdxLxSq[0])*rdxUy[1]+(8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxLySq[1]+(8646.397631383834*rdxLxSq[0]-8646.397631383834*rdxUxSq[0])*rdxLy[1]-360.2665679743264*rdxUxCu[0]-23777.59348630554*rdxLx[0]*rdxUxSq[0]+23777.59348630554*rdxLxSq[0]*rdxUx[0]+360.2665679743264*rdxLxCu[0])*rho[2]+((-360.2665679743264*rdxUyCu[1])+((-23777.59348630554*rdxLy[1])-8646.397631383834*rdxUx[0]-8646.397631383834*rdxLx[0])*rdxUySq[1]+(23777.59348630554*rdxLySq[1]-8286.131063409508*rdxUxSq[0]+6845.064791512203*rdxLx[0]*rdxUx[0]-8286.131063409508*rdxLxSq[0])*rdxUy[1]+360.2665679743264*rdxLyCu[1]+(8646.397631383834*rdxUx[0]+8646.397631383834*rdxLx[0])*rdxLySq[1]+(8286.131063409508*rdxUxSq[0]-6845.064791512203*rdxLx[0]*rdxUx[0]+8286.131063409508*rdxLxSq[0])*rdxLy[1])*rho[1]+(21632.0*rdxUx[0]-21632.0*rdxLx[0])*rho[0]*rdxUySq[1]+(21632.0*rdxUxSq[0]-21632.0*rdxLxSq[0])*rho[0]*rdxUy[1]+(21632.0*rdxLx[0]-21632.0*rdxUx[0])*rho[0]*rdxLySq[1]+(21632.0*rdxLxSq[0]-21632.0*rdxUxSq[0])*rho[0]*rdxLy[1])*volFac+(132.0*rdxUyR4[1]+(8586.0*rdxLy[1]+3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxUyCu[1]+((-17154.0*rdxLySq[1])+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLy[1]+2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-258.0*rdxLyCu[1])+((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1]-63.0*rdxUxCu[0]-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-63.0*rdxUx[0]*rdxUyCu[1])+((-4284.0*rdxUx[0]*rdxLy[1])+2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-4284.0*rdxUx[0]*rdxLySq[1])+((-17516.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0])*rdxLy[1]+3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-63.0*rdxUx[0]*rdxLyCu[1]+(2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+132.0*rdxUxR4[0]+8586.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]-258.0*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((-258.0*rdxLy[1]*rdxUyCu[1])+(((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLy[1]-17154.0*rdxLySq[1])*rdxUySq[1]+(8586.0*rdxLyCu[1]+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+132.0*rdxLyR4[1]+(3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxLyCu[1]+(2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxLySq[1]+((-63.0*rdxUxCu[0])-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-63.0*rdxLx[0]*rdxUyCu[1])+((-4284.0*rdxLx[0]*rdxLy[1])-6158.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-4284.0*rdxLx[0]*rdxLySq[1])+((-10106.0*rdxLx[0]*rdxUx[0])-17516.0*rdxLxSq[0])*rdxLy[1]-6353.0*rdxLx[0]*rdxUxSq[0]-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxUy[1]-63.0*rdxLx[0]*rdxLyCu[1]+(2812.0*rdxLxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6353.0*rdxLx[0]*rdxUxSq[0])-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxLy[1]-258.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]+8586.0*rdxLxCu[0]*rdxUx[0]+132.0*rdxLxR4[0])*phiLx[3]+((8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxUyCu[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLy[1]-7925.864495435182*rdxUxSq[0]+7925.864495435182*rdxLxSq[0])*rdxUySq[1]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1]+157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(77.94228634059945*rdxUx[0]*rdxUyCu[1]+(5300.075471160763*rdxUx[0]*rdxLy[1]-2591.14800812304*rdxUxSq[0]+6730.749438212657*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(5300.075471160763*rdxUx[0]*rdxLySq[1]+(20937.03016189259*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0])*rdxLy[1]-2793.797952608599*rdxUxCu[0]+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+77.94228634059945*rdxUx[0]*rdxLyCu[1]+(6730.749438212657*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxUxSq[0])*rdxLySq[1]+((-2793.797952608599*rdxUxCu[0])+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-124.7076581449591*rdxUxR4[0]-8074.820864886104*rdxLx[0]*rdxUxCu[0]+18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+280.592230826158*rdxLxCu[0]*rdxUx[0])*phiUx[2]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxLyCu[1]+(7925.864495435182*rdxLxSq[0]-7925.864495435182*rdxUxSq[0])*rdxLySq[1]+(157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-77.94228634059945*rdxLx[0]*rdxUyCu[1])+((-5300.075471160763*rdxLx[0]*rdxLy[1])-6730.749438212657*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*rdxUySq[1]+((-5300.075471160763*rdxLx[0]*rdxLySq[1])+((-13236.33227144136*rdxLx[0]*rdxUx[0])-20937.03016189259*rdxLxSq[0])*rdxLy[1]-6933.399382698215*rdxLx[0]*rdxUxSq[0]-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxUy[1]-77.94228634059945*rdxLx[0]*rdxLyCu[1]+(2591.14800812304*rdxLxSq[0]-6730.749438212657*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6933.399382698215*rdxLx[0]*rdxUxSq[0])-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxLy[1]-280.592230826158*rdxLx[0]*rdxUxCu[0]-18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+8074.820864886104*rdxLxCu[0]*rdxUx[0]+124.7076581449591*rdxLxR4[0])*phiLx[2]-124.7076581449591*phiUy[1]*rdxUyR4[1]+(((-8074.820864886104*phiUy[1])-280.592230826158*phiLy[1])*rdxLy[1]+((-2793.797952608599*rdxUx[0])-2793.797952608599*rdxLx[0])*phiUy[1]+157.6166234887678*rdxUx[0]*phiUx[1]+157.6166234887678*rdxLx[0]*phiLx[1]+(7683.0*phiUy[0]-195.0*phiUx[0])*rdxUx[0]+(195.0*phiLx[0]-7683.0*phiUy[0])*rdxLx[0])*rdxUyCu[1]+((18674.97180720763*phiUy[1]-18674.97180720763*phiLy[1])*rdxLySq[1]+((17086.68121666697*rdxUx[0]+17086.68121666697*rdxLx[0])*phiUy[1]+10402.69715025867*rdxUx[0]*phiUx[1]+((-6933.399382698215*rdxUx[0])-6933.399382698215*rdxLx[0])*phiLy[1]+10402.69715025867*rdxLx[0]*phiLx[1]+(3705.0*phiUy[0]-12870.0*phiUx[0]+16653.0*phiLy[0])*rdxUx[0]+((-3705.0*phiUy[0])-16653.0*phiLy[0]+12870.0*phiLx[0])*rdxLx[0])*rdxLy[1]+((-2591.14800812304*rdxUxSq[0])+20937.03016189259*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxLxSq[0])*phiUy[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxUxSq[0])*phiUx[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxLxSq[0])*phiLx[1]+(7488.0*phiUy[0]+7488.0*phiUx[0])*rdxUxSq[0]+(16848.0*phiLx[0]-16848.0*phiUx[0])*rdxLx[0]*rdxUx[0]+((-7488.0*phiUy[0])-7488.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((280.592230826158*phiUy[1]+8074.820864886104*phiLy[1])*rdxLyCu[1]+((6933.399382698215*rdxUx[0]+6933.399382698215*rdxLx[0])*phiUy[1]-10402.69715025867*rdxUx[0]*phiUx[1]+((-17086.68121666697*rdxUx[0])-17086.68121666697*rdxLx[0])*phiLy[1]-10402.69715025867*rdxLx[0]*phiLx[1]+((-16653.0*phiUy[0])+12870.0*phiUx[0]-3705.0*phiLy[0])*rdxUx[0]+(16653.0*phiUy[0]+3705.0*phiLy[0]-12870.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((6730.749438212657*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0]+6730.749438212657*rdxLxSq[0])*phiUy[1]+((-6730.749438212657*rdxUxSq[0])-13236.33227144136*rdxLx[0]*rdxUx[0]-6730.749438212657*rdxLxSq[0])*phiLy[1]+(16848.0*phiLy[0]-16848.0*phiUy[0])*rdxUxSq[0]+(16848.0*phiUy[0]-16848.0*phiLy[0])*rdxLxSq[0])*rdxLy[1]+(77.94228634059945*rdxUxCu[0]+5300.075471160763*rdxLx[0]*rdxUxSq[0]+5300.075471160763*rdxLxSq[0]*rdxUx[0]+77.94228634059945*rdxLxCu[0])*phiUy[1]+((-8083.48111892395*rdxUxCu[0])-2994.715846286589*rdxLx[0]*rdxUxSq[0]+15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(15333.84579940727*rdxLx[0]*rdxUxSq[0]-2994.715846286589*rdxLxSq[0]*rdxUx[0]-8083.48111892395*rdxLxCu[0])*phiLx[1]+(7683.0*phiUx[0]-195.0*phiUy[0])*rdxUxCu[0]+((-12870.0*phiUy[0])+3705.0*phiUx[0]+16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(12870.0*phiUy[0]-16653.0*phiUx[0]-3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(195.0*phiUy[0]-7683.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+124.7076581449591*phiLy[1]*rdxLyR4[1]+((-157.6166234887678*rdxUx[0]*phiUx[1])+(2793.797952608599*rdxUx[0]+2793.797952608599*rdxLx[0])*phiLy[1]-157.6166234887678*rdxLx[0]*phiLx[1]+(195.0*phiUx[0]-7683.0*phiLy[0])*rdxUx[0]+(7683.0*phiLy[0]-195.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((7925.864495435182*rdxUxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiUx[1]+(2591.14800812304*rdxUxSq[0]-20937.03016189259*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*phiLy[1]+(7925.864495435182*rdxLxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiLx[1]+((-7488.0*phiUx[0])-7488.0*phiLy[0])*rdxUxSq[0]+(16848.0*phiUx[0]-16848.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(7488.0*phiLy[0]+7488.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((8083.48111892395*rdxUxCu[0]+2994.715846286589*rdxLx[0]*rdxUxSq[0]-15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-77.94228634059945*rdxUxCu[0])-5300.075471160763*rdxLx[0]*rdxUxSq[0]-5300.075471160763*rdxLxSq[0]*rdxUx[0]-77.94228634059945*rdxLxCu[0])*phiLy[1]+((-15333.84579940727*rdxLx[0]*rdxUxSq[0])+2994.715846286589*rdxLxSq[0]*rdxUx[0]+8083.48111892395*rdxLxCu[0])*phiLx[1]+(195.0*phiLy[0]-7683.0*phiUx[0])*rdxUxCu[0]+((-3705.0*phiUx[0])+12870.0*phiLy[0]-16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(16653.0*phiUx[0]-12870.0*phiLy[0]+3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(7683.0*phiLx[0]-195.0*phiLy[0])*rdxLxCu[0])*rdxLy[1])/(12.0*rdxUyR4[1]+(1608.0*rdxLy[1]+1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxUyCu[1]+(53892.0*rdxLySq[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLy[1]+2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxUySq[1]+(1608.0*rdxLyCu[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLySq[1]+(69048.0*rdxUxSq[0]+166696.0*rdxLx[0]*rdxUx[0]+69048.0*rdxLxSq[0])*rdxLy[1]+1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxUy[1]+12.0*rdxLyR4[1]+(1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxLyCu[1]+(2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxLySq[1]+(1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxLy[1]+12.0*rdxUxR4[0]+1608.0*rdxLx[0]*rdxUxCu[0]+53892.0*rdxLxSq[0]*rdxUxSq[0]+1608.0*rdxLxCu[0]*rdxUx[0]+12.0*rdxLxR4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[0]*rho[1]+144.0*rho[0]*rdxCp2[1]+1600.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(81.0*phiUy[0]+81.0*phiLy[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+93.53074360871933*rdxCp2[0]*phiUx[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+(900.0*phiUy[0]-54.0*phiUx[0]+900.0*phiLy[0]+864.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiUx[1]+(1020.0*phiUx[0]+3120.0*bcVals[0])*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]+277.1281292110203*rdxCp2[0]*rho[0])*volFac+((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]-360.0*rdxCp2[0]*phiUx[1]+189.0*rdxCp2[0]*phiLy[1]+(155.8845726811989*phiUy[0]+311.7691453623978*phiUx[0]+155.8845726811989*phiLy[0]-1247.076581449591*bcVals[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiUx[1]+(623.5382907247956*phiUx[0]-1247.076581449591*bcVals[0])*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = ((748.2459488697547*rdxCp2[0]*rho[3]+(368.0*rdxCp2[1]+1600.0*rdxCp2[0])*rho[2])*volFac-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(1020.0*rdxCp2Sq[0]-138.0*rdxCp2[0]*rdxCp2[1])*phiUx[2]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(199.1858428704209*phiUy[0]-199.1858428704209*phiLy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiUy[0]-866.0254037844386*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[0]*rho[2])*volFac+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiUx[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+623.5382907247956*rdxCp2Sq[0])*phiUx[2]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiUy[0]-150.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[1]-576.0*rho[0]*rdxCp2[1]-1680.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-324.0*phiUy[0])-324.0*phiLy[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+311.7691453623978*rdxCp2[0]*phiUx[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+((-945.0*phiUy[0])-324.0*phiUx[0]-945.0*phiLy[0]+576.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiUx[1]+(2080.0*bcVals[0]-720.0*phiUx[0])*rdxCp2Sq[0]))/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]-138.5640646055102*rdxCp2[0]*rho[0])*volFac+((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]+75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiUx[1]+54.0*rdxCp2[0]*phiLy[1]+((-77.94228634059945*phiUy[0])+155.8845726811989*phiUx[0]-77.94228634059945*phiLy[0]+277.1281292110203*bcVals[0])*rdxCp2[0])*rdxCp2[1]+277.1281292110203*bcVals[0]*rdxCp2Sq[0])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[3]+((-1472.0*rdxCp2[1])-1680.0*rdxCp2[0])*rho[2])*volFac-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-828.0*rdxCp2[0]*rdxCp2[1])-720.0*rdxCp2Sq[0])*phiUx[2]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(796.7433714816835*phiLy[0]-796.7433714816835*phiUy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiLy[0]-909.3266739736605*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[0]*rho[2])*volFac+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]+1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiUx[2]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiLy[0]-225.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[1]-144.0*rho[0]*rdxCp2[1]-1600.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-81.0*phiUy[0])-81.0*phiLy[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+93.53074360871933*rdxCp2[0]*phiLx[1]-864.0*rdxCp2[0]*bcVals[1]+((-900.0*phiUy[0])-900.0*phiLy[0]+54.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiLx[1]-3120.0*rdxCp2Sq[0]*bcVals[1]-1020.0*phiLx[0]*rdxCp2Sq[0]))/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]-277.1281292110203*rdxCp2[0]*rho[0])*volFac+((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]+189.0*rdxCp2[0]*phiLy[1]-360.0*rdxCp2[0]*phiLx[1]+1247.076581449591*rdxCp2[0]*bcVals[1]+((-155.8845726811989*phiUy[0])-155.8845726811989*phiLy[0]-311.7691453623978*phiLx[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiLx[1]+1247.076581449591*rdxCp2Sq[0]*bcVals[1]-623.5382907247956*phiLx[0]*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[3]+((-368.0*rdxCp2[1])-1600.0*rdxCp2[0])*rho[2])*volFac-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiLx[3]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(138.0*rdxCp2[0]*rdxCp2[1]-1020.0*rdxCp2Sq[0])*phiLx[2]+(199.1858428704209*phiLy[0]-199.1858428704209*phiUy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiLy[0]-866.0254037844386*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[0]*rho[2])*volFac+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiLx[3]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-796.7433714816835*rdxCp2[0]*rdxCp2[1])-623.5382907247956*rdxCp2Sq[0])*phiLx[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiLy[0]-150.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[0]*rho[1]+576.0*rho[0]*rdxCp2[1]+1680.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(324.0*phiUy[0]+324.0*phiLy[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+311.7691453623978*rdxCp2[0]*phiLx[1]+576.0*rdxCp2[0]*bcVals[1]+(945.0*phiUy[0]+945.0*phiLy[0]+324.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiLx[1]+2080.0*rdxCp2Sq[0]*bcVals[1]+720.0*phiLx[0]*rdxCp2Sq[0])/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]+138.5640646055102*rdxCp2[0]*rho[0])*volFac+((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]-75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]+54.0*rdxCp2[0]*phiLy[1]-150.0*rdxCp2[0]*phiLx[1]+277.1281292110203*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiUy[0]+77.94228634059945*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1]+277.1281292110203*rdxCp2Sq[0]*bcVals[1])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = ((277.1281292110203*rdxCp2[0]*rho[3]+(1472.0*rdxCp2[1]+1680.0*rdxCp2[0])*rho[2])*volFac-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiLx[3]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(828.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiLx[2]+(796.7433714816835*phiUy[0]-796.7433714816835*phiLy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiUy[0]-909.3266739736605*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[0]*rho[2])*volFac+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]-1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiUy[0]-225.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[1]*rho[2]+1600.0*rho[0]*rdxCp2[1]+144.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(3120.0*rdxCp2Sq[1]+864.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+1020.0*phiUy[0]*rdxCp2Sq[1]+((-866.0254037844386*rdxCp2[0]*phiUx[1])+866.0254037844386*rdxCp2[0]*phiLx[1]+((-54.0*phiUy[0])+900.0*phiUx[0]+900.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-77.94228634059945*rdxCp2Sq[0]*phiUx[1]+77.94228634059945*rdxCp2Sq[0]*phiLx[1]+(81.0*phiUx[0]+81.0*phiLx[0])*rdxCp2Sq[0])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = ((748.2459488697547*rdxCp2[1]*rho[3]+(1600.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[1])*volFac+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+1020.0*phiUy[1]*rdxCp2Sq[1]+((-138.0*rdxCp2[0]*phiUy[1])-700.0*rdxCp2[0]*phiUx[1]-700.0*rdxCp2[0]*phiLx[1]+(866.0254037844386*phiUx[0]-866.0254037844386*phiLx[0])*rdxCp2[0])*rdxCp2[1]-161.0*rdxCp2Sq[0]*phiUx[1]-161.0*rdxCp2Sq[0]*phiLx[1]+(199.1858428704209*phiUx[0]-199.1858428704209*phiLx[0])*rdxCp2Sq[0])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]+277.1281292110203*rho[0]*rdxCp2[1])*volFac+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]+((-1247.076581449591*rdxCp2Sq[1])-1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[2]+623.5382907247956*phiUy[0]*rdxCp2Sq[1]+((-150.0*rdxCp2[0]*phiUx[1])+150.0*rdxCp2[0]*phiLx[1]+(311.7691453623978*phiUy[0]+155.8845726811989*phiUx[0]+155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[1]*rho[1])*volFac+((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]+623.5382907247956*phiUy[1]*rdxCp2Sq[1]+(796.7433714816835*rdxCp2[0]*phiUy[1]-121.2435565298214*rdxCp2[0]*phiUx[1]-121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiUx[0]-150.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[2]-1680.0*rho[0]*rdxCp2[1]-576.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]-720.0*phiUy[0]*rdxCp2Sq[1]+(909.3266739736605*rdxCp2[0]*phiUx[1]-909.3266739736605*rdxCp2[0]*phiLx[1]+((-324.0*phiUy[0])-945.0*phiUx[0]-945.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+311.7691453623978*rdxCp2Sq[0]*phiUx[1]-311.7691453623978*rdxCp2Sq[0]*phiLx[1]+((-324.0*phiUx[0])-324.0*phiLx[0])*rdxCp2Sq[0]))/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[3]+((-1680.0*rdxCp2[1])-1472.0*rdxCp2[0])*rho[1])*volFac+(692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiUy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-720.0*phiUy[1]*rdxCp2Sq[1]+((-828.0*rdxCp2[0]*phiUy[1])+735.0*rdxCp2[0]*phiUx[1]+735.0*rdxCp2[0]*phiLx[1]+(909.3266739736605*phiLx[0]-909.3266739736605*phiUx[0])*rdxCp2[0])*rdxCp2[1]+644.0*rdxCp2Sq[0]*phiUx[1]+644.0*rdxCp2Sq[0]*phiLx[1]+(796.7433714816835*phiLx[0]-796.7433714816835*phiUx[0])*rdxCp2Sq[0]))/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]-138.5640646055102*rho[0]*rdxCp2[1])*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(75.0*rdxCp2[0]*phiUx[1]-75.0*rdxCp2[0]*phiLx[1]+(155.8845726811989*phiUy[0]-77.94228634059945*phiUx[0]-77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[1]*rho[1])*volFac-1150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+(1195.115057222525*rdxCp2[0]*phiUy[1]+181.8653347947321*rdxCp2[0]*phiUx[1]+181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiLx[0]-225.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[2]-1600.0*rho[0]*rdxCp2[1]-144.0*rdxCp2[0]*rho[0])*volFac-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-3120.0*rdxCp2Sq[1])-864.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]-1020.0*phiLy[0]*rdxCp2Sq[1]+(866.0254037844386*rdxCp2[0]*phiUx[1]-866.0254037844386*rdxCp2[0]*phiLx[1]+((-900.0*phiUx[0])+54.0*phiLy[0]-900.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0]*phiUx[1]-77.94228634059945*rdxCp2Sq[0]*phiLx[1]+((-81.0*phiUx[0])-81.0*phiLx[0])*rdxCp2Sq[0]))/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[3]+((-1600.0*rdxCp2[1])-368.0*rdxCp2[0])*rho[1])*volFac-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-1020.0*phiLy[1]*rdxCp2Sq[1]+(700.0*rdxCp2[0]*phiUx[1]+138.0*rdxCp2[0]*phiLy[1]+700.0*rdxCp2[0]*phiLx[1]+(866.0254037844386*phiLx[0]-866.0254037844386*phiUx[0])*rdxCp2[0])*rdxCp2[1]+161.0*rdxCp2Sq[0]*phiUx[1]+161.0*rdxCp2Sq[0]*phiLx[1]+(199.1858428704209*phiLx[0]-199.1858428704209*phiUx[0])*rdxCp2Sq[0]))/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]-277.1281292110203*rho[0]*rdxCp2[1])*volFac+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+(1247.076581449591*rdxCp2Sq[1]+1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]-623.5382907247956*phiLy[0]*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiLx[1]+((-155.8845726811989*phiUx[0])-311.7691453623978*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[1]*rho[1])*volFac+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]-623.5382907247956*phiLy[1]*rdxCp2Sq[1]+(121.2435565298214*rdxCp2[0]*phiUx[1]-796.7433714816835*rdxCp2[0]*phiLy[1]+121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiLx[0]-150.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[1]*rho[2]+1680.0*rho[0]*rdxCp2[1]+576.0*rdxCp2[0]*rho[0])*volFac-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+720.0*phiLy[0]*rdxCp2Sq[1]+((-909.3266739736605*rdxCp2[0]*phiUx[1])+909.3266739736605*rdxCp2[0]*phiLx[1]+(945.0*phiUx[0]+324.0*phiLy[0]+945.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-311.7691453623978*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLx[1]+(324.0*phiUx[0]+324.0*phiLx[0])*rdxCp2Sq[0])/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = ((277.1281292110203*rdxCp2[1]*rho[3]+(1680.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[1])*volFac-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3]+(692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiLy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+720.0*phiLy[1]*rdxCp2Sq[1]+((-735.0*rdxCp2[0]*phiUx[1])+828.0*rdxCp2[0]*phiLy[1]-735.0*rdxCp2[0]*phiLx[1]+(909.3266739736605*phiUx[0]-909.3266739736605*phiLx[0])*rdxCp2[0])*rdxCp2[1]-644.0*rdxCp2Sq[0]*phiUx[1]-644.0*rdxCp2Sq[0]*phiLx[1]+(796.7433714816835*phiUx[0]-796.7433714816835*phiLx[0])*rdxCp2Sq[0])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]+138.5640646055102*rho[0]*rdxCp2[1])*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+((-75.0*rdxCp2[0]*phiUx[1])+75.0*rdxCp2[0]*phiLx[1]+(77.94228634059945*phiUx[0]-155.8845726811989*phiLy[0]+77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[1]*rho[1])*volFac+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+((-181.8653347947321*rdxCp2[0]*phiUx[1])-1195.115057222525*rdxCp2[0]*phiLy[1]-181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiUx[0]-225.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+328950.0*phiUy[0]*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiUx[1]+(2715915.0*phiUy[0]-193500.0*phiUx[0]+3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiUx[1]+((-193500.0*phiUy[0])+2715915.0*phiUx[0]+9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiUx[1]+(328950.0*phiUx[0]+1006200.0*bcVals[0])*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+328950.0*phiUy[1]*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiUx[1]+(567939.4598018348*phiUy[0]+1117172.770881926*phiUx[0]-4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiUx[1]+((-33515.18312645776*phiUy[0])+1919718.512568964*phiUx[0]-4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiUx[1]+(201091.0987587466*phiUx[0]-402182.1975174932*bcVals[0])*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiUx[1]+(1919718.512568964*phiUy[0]-33515.18312645776*phiUx[0]+536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiUx[1]+(1117172.770881926*phiUy[0]+567939.4598018348*phiUx[0]+1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiUx[1]+(44400.0*phiUy[0]+64500.0*phiUx[0]-258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiUx[1]+(64500.0*phiUy[0]+44400.0*phiUx[0]+40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*volFac+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-12960.0*phiUy[0]*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiUx[1]+((-456408.0*phiUy[0])+11340.0*phiUx[0]-181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiUx[1]+((-1842246.0*phiUy[0])-532953.0*phiUx[0]-2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiUx[1]+((-928800.0*phiUy[0])-2563956.0*phiUx[0]-8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiUx[1]+((-1052640.0*phiUx[0])-3219840.0*bcVals[0])*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*volFac+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-12960.0*phiUy[1]*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiUx[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiUx[1]+((-301792.5327108011*phiUy[0])-659547.6270141526*phiUx[0]+1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiUx[1]+((-160872.8790069972*phiUy[0])-1745283.675738703*phiUx[0]+3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiUx[1]+(1286983.032055978*bcVals[0]-643491.516027989*phiUx[0])*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiUx[1]+(98207.28078915528*phiUy[0]+14029.6115413079*phiUx[0]-224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiUx[1]+(2061183.762277152*phiUy[0]-814886.6036909672*phiUx[0]-2492594.31717237*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiUx[1]+(6703036.625291553*phiUy[0]-3407636.758811008*phiUx[0]-6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiUx[1]+(332100.0*bcVals[0]-166050.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiUx[1]+(387000.0*phiUy[0]-266400.0*phiUx[0]-241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*volFac+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-1052640.0*phiUy[0]*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiUx[1]+((-2563956.0*phiUy[0])-928800.0*phiUx[0]+1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiUx[1]+((-532953.0*phiUy[0])-1842246.0*phiUx[0]+5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiUx[1]+(11340.0*phiUy[0]-456408.0*phiUx[0]+1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiUx[1]+(37440.0*bcVals[0]-12960.0*phiUx[0])*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+5263200.0*phiUy[1]*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiUx[1]+((-3407636.758811008*phiUy[0])+6703036.625291553*phiUx[0]+1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiUx[1]+((-814886.6036909672*phiUy[0])+2061183.762277152*phiUx[0]+1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiUx[1]+(14029.6115413079*phiUy[0]+98207.28078915528*phiUx[0]+2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+74824.59488697546*bcVals[0]*rdxCp2R4[0])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiUx[1]+((-1745283.675738703*phiUy[0])-160872.8790069972*phiUx[0]+285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiUx[1]+((-659547.6270141526*phiUy[0])-301792.5327108011*phiUx[0]+871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiUx[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiUx[1]+((-266400.0*phiUy[0])+387000.0*phiUx[0]+688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+(580800.0*bcVals[0]-166050.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*volFac+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+1440.0*phiUy[0]*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiUx[1]+(11799.0*phiUy[0]+1890.0*phiUx[0]-3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiUx[1]+(1890.0*phiUy[0]+11799.0*phiUx[0]-33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiUx[1]+(1440.0*phiUx[0]-4160.0*bcVals[0])*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-1440.0*phiUy[1]*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiUx[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiUx[0]-4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiUx[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiUx[0]-14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1662.768775266122*bcVals[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiUx[1]+((-467.6537180435967*phiUy[0])+467.6537180435967*phiUx[0]-831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiUx[1]+((-2727.980021920981*phiUy[0])+2727.980021920981*phiUx[0]-8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiUx[1]+(1125.0*phiUy[0]-1125.0*phiUx[0]-2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiUx[1]+((-1125.0*phiUy[0])+1125.0*phiUx[0]-4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-328950.0*phiLy[0]*rdxCp2R3[1]+((-335151.8312645776*rdxCp2[0]*phiUx[1])-1533436.541464953*rdxCp2[0]*phiLy[1]+(193500.0*phiUx[0]-2715915.0*phiLy[0]-3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2176434.423012786*rdxCp2Sq[0]*phiUx[1]+90490.99444143593*rdxCp2Sq[0]*phiLy[1]+((-2715915.0*phiUx[0])+193500.0*phiLy[0]-9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiUx[1]+((-328950.0*phiUx[0])-1006200.0*bcVals[0])*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]-328950.0*phiLy[1]*rdxCp2R3[1]+(1290000.0*rdxCp2[0]*phiUx[1]-125505.0*rdxCp2[0]*phiLy[1]+((-1117172.770881926*phiUx[0])-567939.4598018348*phiLy[0]+4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2054550.0*rdxCp2Sq[0]*phiUx[1]+40635.0*rdxCp2Sq[0]*phiLy[1]+((-1919718.512568964*phiUx[0])+33515.18312645776*phiLy[0]+4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiUx[1]+(402182.1975174932*bcVals[0]-201091.0987587466*phiUx[0])*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-58050.0*rdxCp2[0]*phiUx[1])-359640.0*rdxCp2[0]*phiLy[1]+(33515.18312645776*phiUx[0]-1919718.512568964*phiLy[0]-536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(493650.0*rdxCp2Sq[0]*phiUx[1]-522450.0*rdxCp2Sq[0]*phiLy[1]+((-567939.4598018348*phiUx[0])-1117172.770881926*phiLy[0]-1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+(74478.18472546172*rdxCp2[0]*phiUx[1]-437394.7904353685*rdxCp2[0]*phiLy[1]+((-64500.0*phiUx[0])-44400.0*phiLy[0]+258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(40789.79651824706*rdxCp2Sq[0]*phiUx[1]-78202.0939617348*rdxCp2Sq[0]*phiLy[1]+((-44400.0*phiUx[0])-64500.0*phiLy[0]-40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*volFac+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+12960.0*phiLy[0]*rdxCp2R4[1]+(19641.45615783106*rdxCp2[0]*phiUx[1]+176773.1054204795*rdxCp2[0]*phiLy[1]+((-11340.0*phiUx[0])+456408.0*phiLy[0]+181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-386469.0325912283*rdxCp2Sq[0]*phiUx[1])+814839.8383191626*rdxCp2Sq[0]*phiLy[1]+(532953.0*phiUx[0]+1842246.0*phiLy[0]+2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2064285.865273508*rdxCp2R3[0]*phiUx[1])+434356.7733188925*rdxCp2R3[0]*phiLy[1]+(2563956.0*phiUx[0]+928800.0*phiLy[0]+8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiUx[1]+(1052640.0*phiUx[0]+3219840.0*bcVals[0])*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*volFac+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+12960.0*phiLy[1]*rdxCp2R4[1]+((-75600.0*rdxCp2[0]*phiUx[1])+157788.0*rdxCp2[0]*phiLy[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-727155.0*rdxCp2Sq[0]*phiUx[1])+465750.0*rdxCp2Sq[0]*phiLy[1]+(659547.6270141526*phiUx[0]+301792.5327108011*phiLy[0]-1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1862820.0*rdxCp2R3[0]*phiUx[1])+195048.0*rdxCp2R3[0]*phiLy[1]+(1745283.675738703*phiUx[0]+160872.8790069972*phiLy[0]-3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiUx[1]+(643491.516027989*phiUx[0]-1286983.032055978*bcVals[0])*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(24300.0*rdxCp2[0]*phiUx[1]+218700.0*rdxCp2[0]*phiLy[1]+((-14029.6115413079*phiUx[0])-98207.28078915528*phiLy[0]+224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((814886.6036909672*phiUx[0]-2061183.762277152*phiLy[0]+2492594.31717237*bcVals[0])*rdxCp2Sq[0]-664200.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-2961900.0*rdxCp2R3[0]*phiUx[1])-3134700.0*rdxCp2R3[0]*phiLy[1]+(3407636.758811008*phiUx[0]-6703036.625291553*phiLy[0]+6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-31176.91453623978*rdxCp2[0]*phiUx[1])-155884.5726811989*rdxCp2[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-175759.8556980518*rdxCp2Sq[0]*phiUx[1])-687061.2540923841*rdxCp2Sq[0]*phiLy[1]+(166050.0*phiUx[0]-332100.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-244738.7791094823*rdxCp2R3[0]*phiUx[1])-469212.5637704087*rdxCp2R3[0]*phiLy[1]+(266400.0*phiUx[0]-387000.0*phiLy[0]+241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*volFac+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+1052640.0*phiLy[0]*rdxCp2R4[1]+((-893738.2167055405*rdxCp2[0]*phiUx[1])-454351.5678414679*rdxCp2[0]*phiLy[1]+(928800.0*phiUx[0]+2563956.0*phiLy[0]-1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-1772702.040022519*rdxCp2Sq[0]*phiUx[1])-108651.5471587956*rdxCp2Sq[0]*phiLy[1]+(1842246.0*phiUx[0]+532953.0*phiLy[0]-5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-439178.8027671645*rdxCp2R3[0]*phiUx[1])+1870.614872174387*rdxCp2R3[0]*phiLy[1]+(456408.0*phiUx[0]-11340.0*phiLy[0]-1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiUx[1]+(12960.0*phiUx[0]-37440.0*bcVals[0])*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-5263200.0*phiLy[1]*rdxCp2R4[1]+(6450000.0*rdxCp2[0]*phiUx[1]-1408860.0*rdxCp2[0]*phiLy[1]+((-6703036.625291553*phiUx[0])+3407636.758811008*phiLy[0]-1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(1983375.0*rdxCp2Sq[0]*phiUx[1]+63990.0*rdxCp2Sq[0]*phiLy[1]+((-2061183.762277152*phiUx[0])+814886.6036909672*phiLy[0]-1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(94500.0*rdxCp2R3[0]*phiUx[1]+9720.0*rdxCp2R3[0]*phiLy[1]+((-98207.28078915528*phiUx[0])-14029.6115413079*phiLy[0]-2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-74824.59488697546*bcVals[0]*rdxCp2R4[0]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+643491.516027989*phiLy[0]*rdxCp2R4[1]+((-154800.0*rdxCp2[0]*phiUx[1])-106560.0*rdxCp2[0]*phiLy[1]+(160872.8790069972*phiUx[0]+1745283.675738703*phiLy[0]-285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-290400.0*rdxCp2Sq[0]*phiUx[1])-66420.0*rdxCp2Sq[0]*phiLy[1]+(301792.5327108011*phiUx[0]+659547.6270141526*phiLy[0]-871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-63000.0*rdxCp2R3[0]*phiUx[1])-10800.0*rdxCp2R3[0]*phiLy[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+(372390.9236273086*rdxCp2[0]*phiUx[1]-2016730.678300898*rdxCp2[0]*phiLy[1]+((-387000.0*phiUx[0])+266400.0*phiLy[0]-688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((166050.0*phiLy[0]-580800.0*bcVals[0])*rdxCp2Sq[0]-543205.7742697513*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-25980.76211353316*rdxCp2R3[0]*phiUx[1])-18706.14872174387*rdxCp2R3[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*volFac+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-1440.0*phiLy[0]*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiUx[1]+1818.653347947321*rdxCp2[0]*phiLy[1]+((-1890.0*phiUx[0])-11799.0*phiLy[0]+3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(11353.59304361399*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLy[1]+((-11799.0*phiUx[0])-1890.0*phiLy[0]+33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiUx[1]+(4160.0*bcVals[0]-1440.0*phiUx[0])*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+1440.0*phiLy[1]*rdxCp2R3[1]+((-2625.0*rdxCp2[0]*phiUx[1])+2664.0*rdxCp2[0]*phiLy[1]+(2727.980021920981*phiUx[0]-2727.980021920981*phiLy[0]+4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiUx[1])+324.0*rdxCp2Sq[0]*phiLy[1]+(467.6537180435967*phiUx[0]-467.6537180435967*phiLy[0]+14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1662.768775266122*bcVals[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(450.0*rdxCp2[0]*phiUx[1]+450.0*rdxCp2[0]*phiLy[1]+((-467.6537180435967*phiUx[0])+467.6537180435967*phiLy[0]+831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2625.0*rdxCp2Sq[0]*phiUx[1]-450.0*rdxCp2Sq[0]*phiLy[1]+((-2727.980021920981*phiUx[0])+2727.980021920981*phiLy[0]+8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-1082.531754730548*rdxCp2[0]*phiUx[1])-4546.633369868302*rdxCp2[0]*phiLy[1]+(1125.0*phiUx[0]-1125.0*phiLy[0]+2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1082.531754730548*rdxCp2Sq[0]*phiUx[1]-779.4228634059946*rdxCp2Sq[0]*phiLy[1]+((-1125.0*phiUx[0])+1125.0*phiLy[0]+4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-328950.0*phiUy[0]*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiLx[1]-3096000.0*rdxCp2[0]*bcVals[1]+(193500.0*phiLx[0]-2715915.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiLx[1]-9081960.0*rdxCp2Sq[0]*bcVals[1]+(193500.0*phiUy[0]-2715915.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiLx[1]-1006200.0*rdxCp2R3[0]*bcVals[1]-328950.0*phiLx[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+328950.0*phiUy[1]*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiLx[1]+4468691.083527703*rdxCp2[0]*bcVals[1]+((-567939.4598018348*phiUy[0])-1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiLx[1]+4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(33515.18312645776*phiUy[0]-1919718.512568964*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiLx[1]+402182.1975174932*rdxCp2R3[0]*bcVals[1]-201091.0987587466*phiLx[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiLx[1]-536242.9300233242*rdxCp2[0]*bcVals[1]+(33515.18312645776*phiLx[0]-1919718.512568964*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiLx[1]-1098466.622160182*rdxCp2Sq[0]*bcVals[1]+((-1117172.770881926*phiUy[0])-567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiLx[1]+258000.0*rdxCp2[0]*bcVals[1]+((-44400.0*phiUy[0])-64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiLx[1]-40200.0*rdxCp2Sq[0]*bcVals[1]+((-64500.0*phiUy[0])-44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*volFac+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+12960.0*phiUy[0]*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiLx[1]+181440.0*rdxCp2[0]*bcVals[1]+(456408.0*phiUy[0]-11340.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiLx[1]+2626452.0*rdxCp2Sq[0]*bcVals[1]+(1842246.0*phiUy[0]+532953.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiLx[1]+8373744.0*rdxCp2R3[0]*bcVals[1]+(928800.0*phiUy[0]+2563956.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiLx[1]+3219840.0*rdxCp2R4[0]*bcVals[1]+1052640.0*phiLx[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*volFac+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-12960.0*phiUy[1]*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiLx[1]-261886.0821044141*rdxCp2[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiLx[1]-1922680.319449907*rdxCp2Sq[0]*bcVals[1]+(301792.5327108011*phiUy[0]+659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiLx[1]-3812313.1094914*rdxCp2R3[0]*bcVals[1]+(160872.8790069972*phiUy[0]+1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiLx[1]-1286983.032055978*rdxCp2R4[0]*bcVals[1]+643491.516027989*phiLx[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiLx[1]+224473.7846609264*rdxCp2[0]*bcVals[1]+((-98207.28078915528*phiUy[0])-14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiLx[1]+2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiLx[0]-2061183.762277152*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiLx[1]+6590799.732961089*rdxCp2R3[0]*bcVals[1]+(3407636.758811008*phiLx[0]-6703036.625291553*phiUy[0])*rdxCp2R3[0])*rdxCp2[1]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiLx[1]-108000.0*rdxCp2[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiLx[1]-332100.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiLx[1]+241200.0*rdxCp2R3[0]*bcVals[1]+(266400.0*phiLx[0]-387000.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*volFac+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1052640.0*phiUy[0]*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiLx[1]+1651200.0*rdxCp2[0]*bcVals[1]+(2563956.0*phiUy[0]+928800.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiLx[1]+5004704.0*rdxCp2Sq[0]*bcVals[1]+(532953.0*phiUy[0]+1842246.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiLx[1]+1303392.0*rdxCp2R3[0]*bcVals[1]+(456408.0*phiLx[0]-11340.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiLx[1]+37440.0*rdxCp2R4[0]*bcVals[1]+12960.0*phiLx[0]*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+5263200.0*phiUy[1]*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiLx[1]+1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(3407636.758811008*phiUy[0]-6703036.625291553*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiLx[1]+1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiUy[0]-2061183.762277152*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiLx[1]+2731097.713374604*rdxCp2R3[0]*bcVals[1]+((-14029.6115413079*phiUy[0])-98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+74824.59488697546*rdxCp2R4[0]*bcVals[1])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiLx[1]+285996.229345773*rdxCp2[0]*bcVals[1]+(1745283.675738703*phiUy[0]+160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiLx[1]+871845.0944978701*rdxCp2Sq[0]*bcVals[1]+(659547.6270141526*phiUy[0]+301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiLx[1]+201610.7140010173*rdxCp2R3[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiLx[1]+688000.0*rdxCp2[0]*bcVals[1]+(266400.0*phiUy[0]-387000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+580800.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiLx[1]+99600.0*rdxCp2R3[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*volFac+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-1440.0*phiUy[0]*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiLx[1]-3360.0*rdxCp2[0]*bcVals[1]+((-11799.0*phiUy[0])-1890.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiLx[1]-33726.0*rdxCp2Sq[0]*bcVals[1]+((-1890.0*phiUy[0])-11799.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiLx[1]-4160.0*rdxCp2R3[0]*bcVals[1]-1440.0*phiLx[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-1440.0*phiUy[1]*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiLx[1]-4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiLx[1]-14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]-1662.768775266122*rdxCp2R3[0]*bcVals[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiLx[1]-831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiLx[1]-8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiLx[1]-2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiLx[1]-4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiUy[0]-1125.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*volFac+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+328950.0*phiLy[0]*rdxCp2R3[1]+((-1533436.541464953*rdxCp2[0]*phiLy[1])-335151.8312645776*rdxCp2[0]*phiLx[1]+3096000.0*rdxCp2[0]*bcVals[1]+(2715915.0*phiLy[0]-193500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(90490.99444143593*rdxCp2Sq[0]*phiLy[1]+2176434.423012786*rdxCp2Sq[0]*phiLx[1]+9081960.0*rdxCp2Sq[0]*bcVals[1]+(2715915.0*phiLx[0]-193500.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiLx[1]+1006200.0*rdxCp2R3[0]*bcVals[1]+328950.0*phiLx[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*volFac+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]-328950.0*phiLy[1]*rdxCp2R3[1]+((-125505.0*rdxCp2[0]*phiLy[1])+1290000.0*rdxCp2[0]*phiLx[1]-4468691.083527703*rdxCp2[0]*bcVals[1]+(567939.4598018348*phiLy[0]+1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(40635.0*rdxCp2Sq[0]*phiLy[1]+2054550.0*rdxCp2Sq[0]*phiLx[1]-4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(1919718.512568964*phiLx[0]-33515.18312645776*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiLx[1]-402182.1975174932*rdxCp2R3[0]*bcVals[1]+201091.0987587466*phiLx[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-359640.0*rdxCp2[0]*phiLy[1])-58050.0*rdxCp2[0]*phiLx[1]+536242.9300233242*rdxCp2[0]*bcVals[1]+(1919718.512568964*phiLy[0]-33515.18312645776*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-522450.0*rdxCp2Sq[0]*phiLy[1])+493650.0*rdxCp2Sq[0]*phiLx[1]+1098466.622160182*rdxCp2Sq[0]*bcVals[1]+(1117172.770881926*phiLy[0]+567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+((-437394.7904353685*rdxCp2[0]*phiLy[1])+74478.18472546172*rdxCp2[0]*phiLx[1]-258000.0*rdxCp2[0]*bcVals[1]+(44400.0*phiLy[0]+64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-78202.0939617348*rdxCp2Sq[0]*phiLy[1])+40789.79651824706*rdxCp2Sq[0]*phiLx[1]+40200.0*rdxCp2Sq[0]*bcVals[1]+(64500.0*phiLy[0]+44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*volFac+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-12960.0*phiLy[0]*rdxCp2R4[1]+(176773.1054204795*rdxCp2[0]*phiLy[1]+19641.45615783106*rdxCp2[0]*phiLx[1]-181440.0*rdxCp2[0]*bcVals[1]+(11340.0*phiLx[0]-456408.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+(814839.8383191626*rdxCp2Sq[0]*phiLy[1]-386469.0325912283*rdxCp2Sq[0]*phiLx[1]-2626452.0*rdxCp2Sq[0]*bcVals[1]+((-1842246.0*phiLy[0])-532953.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(434356.7733188925*rdxCp2R3[0]*phiLy[1]-2064285.865273508*rdxCp2R3[0]*phiLx[1]-8373744.0*rdxCp2R3[0]*bcVals[1]+((-928800.0*phiLy[0])-2563956.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiLx[1]-3219840.0*rdxCp2R4[0]*bcVals[1]-1052640.0*phiLx[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*volFac+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+12960.0*phiLy[1]*rdxCp2R4[1]+(157788.0*rdxCp2[0]*phiLy[1]-75600.0*rdxCp2[0]*phiLx[1]+261886.0821044141*rdxCp2[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(465750.0*rdxCp2Sq[0]*phiLy[1]-727155.0*rdxCp2Sq[0]*phiLx[1]+1922680.319449907*rdxCp2Sq[0]*bcVals[1]+((-301792.5327108011*phiLy[0])-659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(195048.0*rdxCp2R3[0]*phiLy[1]-1862820.0*rdxCp2R3[0]*phiLx[1]+3812313.1094914*rdxCp2R3[0]*bcVals[1]+((-160872.8790069972*phiLy[0])-1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiLx[1]+1286983.032055978*rdxCp2R4[0]*bcVals[1]-643491.516027989*phiLx[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+(218700.0*rdxCp2[0]*phiLy[1]+24300.0*rdxCp2[0]*phiLx[1]-224473.7846609264*rdxCp2[0]*bcVals[1]+(98207.28078915528*phiLy[0]+14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-664200.0*rdxCp2Sq[0]*phiLx[1])-2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLy[0]-814886.6036909672*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3134700.0*rdxCp2R3[0]*phiLy[1])-2961900.0*rdxCp2R3[0]*phiLx[1]-6590799.732961089*rdxCp2R3[0]*bcVals[1]+(6703036.625291553*phiLy[0]-3407636.758811008*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-155884.5726811989*rdxCp2[0]*phiLy[1])-31176.91453623978*rdxCp2[0]*phiLx[1]+108000.0*rdxCp2[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-687061.2540923841*rdxCp2Sq[0]*phiLy[1])-175759.8556980518*rdxCp2Sq[0]*phiLx[1]+332100.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-469212.5637704087*rdxCp2R3[0]*phiLy[1])-244738.7791094823*rdxCp2R3[0]*phiLx[1]-241200.0*rdxCp2R3[0]*bcVals[1]+(387000.0*phiLy[0]-266400.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*volFac+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1052640.0*phiLy[0]*rdxCp2R4[1]+((-454351.5678414679*rdxCp2[0]*phiLy[1])-893738.2167055405*rdxCp2[0]*phiLx[1]-1651200.0*rdxCp2[0]*bcVals[1]+((-2563956.0*phiLy[0])-928800.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-108651.5471587956*rdxCp2Sq[0]*phiLy[1])-1772702.040022519*rdxCp2Sq[0]*phiLx[1]-5004704.0*rdxCp2Sq[0]*bcVals[1]+((-532953.0*phiLy[0])-1842246.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(1870.614872174387*rdxCp2R3[0]*phiLy[1]-439178.8027671645*rdxCp2R3[0]*phiLx[1]-1303392.0*rdxCp2R3[0]*bcVals[1]+(11340.0*phiLy[0]-456408.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiLx[1]-37440.0*rdxCp2R4[0]*bcVals[1]-12960.0*phiLx[0]*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*volFac+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]-5263200.0*phiLy[1]*rdxCp2R4[1]+((-1408860.0*rdxCp2[0]*phiLy[1])+6450000.0*rdxCp2[0]*phiLx[1]-1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(6703036.625291553*phiLx[0]-3407636.758811008*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+(63990.0*rdxCp2Sq[0]*phiLy[1]+1983375.0*rdxCp2Sq[0]*phiLx[1]-1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLx[0]-814886.6036909672*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(9720.0*rdxCp2R3[0]*phiLy[1]+94500.0*rdxCp2R3[0]*phiLx[1]-2731097.713374604*rdxCp2R3[0]*bcVals[1]+(14029.6115413079*phiLy[0]+98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-74824.59488697546*rdxCp2R4[0]*bcVals[1]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]-643491.516027989*phiLy[0]*rdxCp2R4[1]+((-106560.0*rdxCp2[0]*phiLy[1])-154800.0*rdxCp2[0]*phiLx[1]-285996.229345773*rdxCp2[0]*bcVals[1]+((-1745283.675738703*phiLy[0])-160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-66420.0*rdxCp2Sq[0]*phiLy[1])-290400.0*rdxCp2Sq[0]*phiLx[1]-871845.0944978701*rdxCp2Sq[0]*bcVals[1]+((-659547.6270141526*phiLy[0])-301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-10800.0*rdxCp2R3[0]*phiLy[1])-63000.0*rdxCp2R3[0]*phiLx[1]-201610.7140010173*rdxCp2R3[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+((-2016730.678300898*rdxCp2[0]*phiLy[1])+372390.9236273086*rdxCp2[0]*phiLx[1]-688000.0*rdxCp2[0]*bcVals[1]+(387000.0*phiLx[0]-266400.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-543205.7742697513*rdxCp2Sq[0]*phiLy[1])-580800.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-18706.14872174387*rdxCp2R3[0]*phiLy[1])-25980.76211353316*rdxCp2R3[0]*phiLx[1]-99600.0*rdxCp2R3[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*volFac+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+1440.0*phiLy[0]*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiLy[1]+1818.653347947321*rdxCp2[0]*phiLx[1]+3360.0*rdxCp2[0]*bcVals[1]+(11799.0*phiLy[0]+1890.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2Sq[0]*phiLy[1]+11353.59304361399*rdxCp2Sq[0]*phiLx[1]+33726.0*rdxCp2Sq[0]*bcVals[1]+(1890.0*phiLy[0]+11799.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiLx[1]+4160.0*rdxCp2R3[0]*bcVals[1]+1440.0*phiLx[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*volFac+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+1440.0*phiLy[1]*rdxCp2R3[1]+(2664.0*rdxCp2[0]*phiLy[1]-2625.0*rdxCp2[0]*phiLx[1]+4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLy[0]-2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(324.0*rdxCp2Sq[0]*phiLy[1]-450.0*rdxCp2Sq[0]*phiLx[1]+14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLy[0]-467.6537180435967*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+1662.768775266122*rdxCp2R3[0]*bcVals[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+(450.0*rdxCp2[0]*phiLy[1]+450.0*rdxCp2[0]*phiLx[1]+831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiLy[1])+2625.0*rdxCp2Sq[0]*phiLx[1]+8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-4546.633369868302*rdxCp2[0]*phiLy[1])-1082.531754730548*rdxCp2[0]*phiLx[1]+2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLy[0]-1125.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-779.4228634059946*rdxCp2Sq[0]*phiLy[1])+1082.531754730548*rdxCp2Sq[0]*phiLx[1]+4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 
  double *dxUy = dx[3]; 
  double *dxLy = dx[4]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdxCp2[2]; 
  double rdxLx[2]; 
  double rdxUx[2]; 
  double rdxLxSq[2]; 
  double rdxUxSq[2]; 
  double rdxLxCu[2]; 
  double rdxUxCu[2]; 
  double rdxLxR4[2]; 
  double rdxUxR4[2]; 
  double rdxLy[2]; 
  double rdxUy[2]; 
  double rdxLySq[2]; 
  double rdxUySq[2]; 
  double rdxLyCu[2]; 
  double rdxUyCu[2]; 
  double rdxLyR4[2]; 
  double rdxUyR4[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = volFac*4.0/(dxLx[1]*dxLx[1]); 
  rdxUx[1]   = volFac*4.0/(dxUx[1]*dxUx[1]); 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = volFac*4.0/(dxLy[0]*dxLy[0]); 
  rdxUy[0]   = volFac*4.0/(dxUy[0]*dxUy[0]); 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = volFac*4.0/(dxLy[1]*dxLy[1]); 
  rdxUy[1]   = volFac*4.0/(dxUy[1]*dxUy[1]); 
  rdxLySq[1] = rdxLy[1]*rdxLy[1]; 
  rdxUySq[1] = rdxUy[1]*rdxUy[1]; 
  rdxLyCu[1] = rdxLySq[1]*rdxLy[1]; 
  rdxUyCu[1] = rdxUySq[1]*rdxUy[1]; 
  rdxLyR4[1] = rdxLyCu[1]*rdxLy[1]; 
  rdxUyR4[1] = rdxUyCu[1]*rdxUy[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((((9600.0*rdxUx[0]-9600.0*rdxLx[0])*rdxUySq[1]+(9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUy[1]+(9600.0*rdxLx[0]-9600.0*rdxUx[0])*rdxLySq[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLy[1])*rho[3]+((-415.6921938165305*rdxUyCu[1])+((-27435.68479189101*rdxLy[1])-25495.78788741387*rdxUx[0]-25495.78788741387*rdxLx[0])*rdxUySq[1]+(27435.68479189101*rdxLySq[1]-25080.09569359734*rdxUxSq[0]-23140.1987891202*rdxLx[0]*rdxUx[0]-25080.09569359734*rdxLxSq[0])*rdxUy[1]+415.6921938165305*rdxLyCu[1]+(25495.78788741387*rdxUx[0]+25495.78788741387*rdxLx[0])*rdxLySq[1]+(25080.09569359734*rdxUxSq[0]+23140.1987891202*rdxLx[0]*rdxUx[0]+25080.09569359734*rdxLxSq[0])*rdxLy[1])*rho[2]+((25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxUySq[1]+((23140.1987891202*rdxLx[0]-23140.1987891202*rdxUx[0])*rdxLy[1]-25495.78788741387*rdxUxSq[0]+25495.78788741387*rdxLxSq[0])*rdxUy[1]+(25080.09569359734*rdxLx[0]-25080.09569359734*rdxUx[0])*rdxLySq[1]+(25495.78788741387*rdxLxSq[0]-25495.78788741387*rdxUxSq[0])*rdxLy[1]-415.6921938165305*rdxUxCu[0]-27435.68479189101*rdxLx[0]*rdxUxSq[0]+27435.68479189101*rdxLxSq[0]*rdxUx[0]+415.6921938165305*rdxLxCu[0])*rho[1]+1104.0*rho[0]*rdxUyCu[1]+(75072.0*rho[0]*rdxLy[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0])*rdxUySq[1]+(75072.0*rho[0]*rdxLySq[1]+(164368.0*rdxUx[0]+164368.0*rdxLx[0])*rho[0]*rdxLy[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0])*rdxUy[1]+1104.0*rho[0]*rdxLyCu[1]+(68144.0*rdxUx[0]+68144.0*rdxLx[0])*rho[0]*rdxLySq[1]+(68144.0*rdxUxSq[0]+164368.0*rdxLx[0]*rdxUx[0]+68144.0*rdxLxSq[0])*rho[0]*rdxLy[1]+(1104.0*rdxUxCu[0]+75072.0*rdxLx[0]*rdxUxSq[0]+75072.0*rdxLxSq[0]*rdxUx[0]+1104.0*rdxLxCu[0])*rho[0])*omega*volFac+(((9375.0*rdxUx[0]-9375.0*rdxLx[0])*rdxUyCu[1]+((12525.0*rdxUx[0]-12525.0*rdxLx[0])*rdxLy[1]+9600.0*rdxUxSq[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+((17775.0*rdxUx[0]-17775.0*rdxLx[0])*rdxLySq[1]+(18000.0*rdxUxSq[0]-18000.0*rdxLxSq[0])*rdxLy[1]+225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]+9600.0*rdxUxSq[0]+18000.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])+9375.0*rdxUxCu[0]+12525.0*rdxLx[0]*rdxUxSq[0]+17775.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+((-9600.0*rdxUxSq[0])-18000.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-9375.0*rdxUxCu[0])-12525.0*rdxLx[0]*rdxUxSq[0]-17775.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((17775.0*rdxLx[0]-17775.0*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((12525.0*rdxLx[0]-12525.0*rdxUx[0])*rdxLySq[1]+(18000.0*rdxLxSq[0]-18000.0*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9375.0*rdxLx[0]-9375.0*rdxUx[0])*rdxLyCu[1]+(9600.0*rdxLxSq[0]-9600.0*rdxUxSq[0])*rdxLySq[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-18000.0*rdxLx[0]*rdxUx[0]-9600.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-17775.0*rdxLx[0]*rdxUxSq[0]-12525.0*rdxLxSq[0]*rdxUx[0]-9375.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(18000.0*rdxLx[0]*rdxUx[0]+9600.0*rdxLxSq[0])*rdxLySq[1]+(17775.0*rdxLx[0]*rdxUxSq[0]+12525.0*rdxLxSq[0]*rdxUx[0]+9375.0*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((-415.6921938165305*rdxUyR4[1])+((-28630.79984911354*rdxLy[1])-25729.61474643567*rdxUx[0]-25729.61474643567*rdxLx[0])*rdxUyCu[1]+((-52637.02404201817*rdxLySq[1])+((-88966.78973077537*rdxUx[0])-88966.78973077537*rdxLx[0])*rdxLy[1]-25911.4800812304*rdxUxSq[0]-78842.95276053529*rdxLx[0]*rdxUx[0]-25911.4800812304*rdxLxSq[0])*rdxUySq[1]+((-779.4228634059946*rdxLyCu[1])+((-48038.4291479228*rdxUx[0])-48038.4291479228*rdxLx[0])*rdxLySq[1]+((-47856.56381312806*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0]-47856.56381312806*rdxLxSq[0])*rdxLy[1]-597.5575286112626*rdxUxCu[0]-40633.91194556586*rdxLx[0]*rdxUxSq[0]-40633.91194556586*rdxLxSq[0]*rdxUx[0]-597.5575286112626*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-233.8268590217983*rdxUx[0]*rdxUyCu[1])+((-15432.57269543869*rdxUx[0]*rdxLy[1])-9145.22826396367*rdxUxSq[0]-19537.53310937693*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15432.57269543869*rdxUx[0]*rdxLySq[1]-8911.401404941873*rdxUxCu[0]-13016.36181888011*rdxLx[0]*rdxUxSq[0]-19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]+19537.53310937693*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8911.401404941873*rdxUxCu[0]+13016.36181888011*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+(779.4228634059946*rdxLy[1]*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(48038.4291479228*rdxUx[0]+48038.4291479228*rdxLx[0])*rdxLy[1])*rdxUySq[1]+(28630.79984911354*rdxLyCu[1]+(88966.78973077537*rdxUx[0]+88966.78973077537*rdxLx[0])*rdxLySq[1]+(47856.56381312806*rdxUxSq[0]+99090.62670101546*rdxLx[0]*rdxUx[0]+47856.56381312806*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+415.6921938165305*rdxLyR4[1]+(25729.61474643567*rdxUx[0]+25729.61474643567*rdxLx[0])*rdxLyCu[1]+(25911.4800812304*rdxUxSq[0]+78842.95276053529*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*rdxLySq[1]+(597.5575286112626*rdxUxCu[0]+40633.91194556586*rdxLx[0]*rdxUxSq[0]+40633.91194556586*rdxLxSq[0]*rdxUx[0]+597.5575286112626*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15432.57269543869*rdxLx[0]*rdxLy[1])-19537.53310937693*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(15432.57269543869*rdxLx[0]*rdxLySq[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-13016.36181888011*rdxLxSq[0]*rdxUx[0]-8911.401404941873*rdxLxCu[0])*rdxUy[1]+233.8268590217983*rdxLx[0]*rdxLyCu[1]+(19537.53310937693*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(19303.70625035514*rdxLx[0]*rdxUxSq[0]+13016.36181888011*rdxLxSq[0]*rdxUx[0]+8911.401404941873*rdxLxCu[0])*rdxLy[1])*phiLx[2]+(396.0*phiUy[0]-36.0*phiC[0])*rdxUyR4[1]+((27378.0*phiUy[0]+846.0*phiLy[0]-4824.0*phiC[0])*rdxLy[1]+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiUy[1]-597.5575286112626*rdxUx[0]*phiUx[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(24531.0*phiUy[0]+621.0*phiUx[0]-3072.0*phiC[0])*rdxUx[0]+(24531.0*phiUy[0]+621.0*phiLx[0]-3072.0*phiC[0])*rdxLx[0])*rdxUyCu[1]+((57078.0*phiUy[0]+57078.0*phiLy[0]-161676.0*phiC[0])*rdxLySq[1]+((13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(92457.0*phiUy[0]+42228.0*phiUx[0]+52131.0*phiLy[0]-208896.0*phiC[0])*rdxUx[0]+(92457.0*phiUy[0]+52131.0*phiLy[0]+42228.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0])*rdxLy[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiUy[1]+((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUy[0]+24756.0*phiUx[0]-6072.0*phiC[0])*rdxUxSq[0]+(79932.0*phiUy[0]+51906.0*phiUx[0]+51906.0*phiLx[0]-207144.0*phiC[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiUy[0]+24756.0*phiLx[0]-6072.0*phiC[0])*rdxLxSq[0])*rdxUySq[1]+((846.0*phiUy[0]+27378.0*phiLy[0]-4824.0*phiC[0])*rdxLyCu[1]+((19303.70625035514*rdxLx[0]-19303.70625035514*rdxUx[0])*phiUy[1]-40633.91194556586*rdxUx[0]*phiUx[1]+(13016.36181888011*rdxLx[0]-13016.36181888011*rdxUx[0])*phiLy[1]+40633.91194556586*rdxLx[0]*phiLx[1]+(52131.0*phiUy[0]+42228.0*phiUx[0]+92457.0*phiLy[0]-208896.0*phiC[0])*rdxUx[0]+(52131.0*phiUy[0]+92457.0*phiLy[0]+42228.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0])*rdxLySq[1]+((19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiUy[1]+((-78842.95276053529*rdxUxSq[0])-99090.62670101546*rdxLx[0]*rdxUx[0])*phiUx[1]+(19537.53310937693*rdxLxSq[0]-19537.53310937693*rdxUxSq[0])*phiLy[1]+(99090.62670101546*rdxLx[0]*rdxUx[0]+78842.95276053529*rdxLxSq[0])*phiLx[1]+(51906.0*phiUy[0]+79932.0*phiUx[0]+51906.0*phiLy[0]-207144.0*phiC[0])*rdxUxSq[0]+(104982.0*phiUy[0]+104982.0*phiUx[0]+104982.0*phiLy[0]+104982.0*phiLx[0]-500088.0*phiC[0])*rdxLx[0]*rdxUx[0]+(51906.0*phiUy[0]+51906.0*phiLy[0]+79932.0*phiLx[0]-207144.0*phiC[0])*rdxLxSq[0])*rdxLy[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiUy[1]+((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(621.0*phiUy[0]+24531.0*phiUx[0]-3072.0*phiC[0])*rdxUxCu[0]+(42228.0*phiUy[0]+92457.0*phiUx[0]+52131.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0]*rdxUxSq[0]+(42228.0*phiUy[0]+52131.0*phiUx[0]+92457.0*phiLx[0]-208896.0*phiC[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiUy[0]+24531.0*phiLx[0]-3072.0*phiC[0])*rdxLxCu[0])*rdxUy[1]+(396.0*phiLy[0]-36.0*phiC[0])*rdxLyR4[1]+((-597.5575286112626*rdxUx[0]*phiUx[1])+(8911.401404941873*rdxLx[0]-8911.401404941873*rdxUx[0])*phiLy[1]+597.5575286112626*rdxLx[0]*phiLx[1]+(621.0*phiUx[0]+24531.0*phiLy[0]-3072.0*phiC[0])*rdxUx[0]+(24531.0*phiLy[0]+621.0*phiLx[0]-3072.0*phiC[0])*rdxLx[0])*rdxLyCu[1]+(((-25911.4800812304*rdxUxSq[0])-47856.56381312806*rdxLx[0]*rdxUx[0])*phiUx[1]+(9145.22826396367*rdxLxSq[0]-9145.22826396367*rdxUxSq[0])*phiLy[1]+(47856.56381312806*rdxLx[0]*rdxUx[0]+25911.4800812304*rdxLxSq[0])*phiLx[1]+(24756.0*phiUx[0]+24756.0*phiLy[0]-6072.0*phiC[0])*rdxUxSq[0]+(51906.0*phiUx[0]+79932.0*phiLy[0]+51906.0*phiLx[0]-207144.0*phiC[0])*rdxLx[0]*rdxUx[0]+(24756.0*phiLy[0]+24756.0*phiLx[0]-6072.0*phiC[0])*rdxLxSq[0])*rdxLySq[1]+(((-25729.61474643567*rdxUxCu[0])-88966.78973077537*rdxLx[0]*rdxUxSq[0]-48038.4291479228*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-233.8268590217983*rdxUxCu[0])-15432.57269543869*rdxLx[0]*rdxUxSq[0]+15432.57269543869*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*phiLy[1]+(48038.4291479228*rdxLx[0]*rdxUxSq[0]+88966.78973077537*rdxLxSq[0]*rdxUx[0]+25729.61474643567*rdxLxCu[0])*phiLx[1]+(24531.0*phiUx[0]+621.0*phiLy[0]-3072.0*phiC[0])*rdxUxCu[0]+(92457.0*phiUx[0]+42228.0*phiLy[0]+52131.0*phiLx[0]-208896.0*phiC[0])*rdxLx[0]*rdxUxSq[0]+(52131.0*phiUx[0]+42228.0*phiLy[0]+92457.0*phiLx[0]-208896.0*phiC[0])*rdxLxSq[0]*rdxUx[0]+(621.0*phiLy[0]+24531.0*phiLx[0]-3072.0*phiC[0])*rdxLxCu[0])*rdxLy[1]+((-415.6921938165305*rdxUxR4[0])-28630.79984911354*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(779.4228634059946*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+28630.79984911354*rdxLxCu[0]*rdxUx[0]+415.6921938165305*rdxLxR4[0])*phiLx[1]+(396.0*phiUx[0]-36.0*phiC[0])*rdxUxR4[0]+(27378.0*phiUx[0]+846.0*phiLx[0]-4824.0*phiC[0])*rdxLx[0]*rdxUxCu[0]+(57078.0*phiUx[0]+57078.0*phiLx[0]-161676.0*phiC[0])*rdxLxSq[0]*rdxUxSq[0]+(846.0*phiUx[0]+27378.0*phiLx[0]-4824.0*phiC[0])*rdxLxCu[0]*rdxUx[0]+(396.0*phiLx[0]-36.0*phiC[0])*rdxLxR4[0])*omega+36.0*phiC[0]*rdxUyR4[1]+(4824.0*phiC[0]*rdxLy[1]+3072.0*phiC[0]*rdxUx[0]+3072.0*phiC[0]*rdxLx[0])*rdxUyCu[1]+(161676.0*phiC[0]*rdxLySq[1]+(208896.0*phiC[0]*rdxUx[0]+208896.0*phiC[0]*rdxLx[0])*rdxLy[1]+6072.0*phiC[0]*rdxUxSq[0]+207144.0*phiC[0]*rdxLx[0]*rdxUx[0]+6072.0*phiC[0]*rdxLxSq[0])*rdxUySq[1]+(4824.0*phiC[0]*rdxLyCu[1]+(208896.0*phiC[0]*rdxUx[0]+208896.0*phiC[0]*rdxLx[0])*rdxLySq[1]+(207144.0*phiC[0]*rdxUxSq[0]+500088.0*phiC[0]*rdxLx[0]*rdxUx[0]+207144.0*phiC[0]*rdxLxSq[0])*rdxLy[1]+3072.0*phiC[0]*rdxUxCu[0]+208896.0*phiC[0]*rdxLx[0]*rdxUxSq[0]+208896.0*phiC[0]*rdxLxSq[0]*rdxUx[0]+3072.0*phiC[0]*rdxLxCu[0])*rdxUy[1]+36.0*phiC[0]*rdxLyR4[1]+(3072.0*phiC[0]*rdxUx[0]+3072.0*phiC[0]*rdxLx[0])*rdxLyCu[1]+(6072.0*phiC[0]*rdxUxSq[0]+207144.0*phiC[0]*rdxLx[0]*rdxUx[0]+6072.0*phiC[0]*rdxLxSq[0])*rdxLySq[1]+(3072.0*phiC[0]*rdxUxCu[0]+208896.0*phiC[0]*rdxLx[0]*rdxUxSq[0]+208896.0*phiC[0]*rdxLxSq[0]*rdxUx[0]+3072.0*phiC[0]*rdxLxCu[0])*rdxLy[1]+36.0*phiC[0]*rdxUxR4[0]+4824.0*phiC[0]*rdxLx[0]*rdxUxCu[0]+161676.0*phiC[0]*rdxLxSq[0]*rdxUxSq[0]+4824.0*phiC[0]*rdxLxCu[0]*rdxUx[0]+36.0*phiC[0]*rdxLxR4[0])/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxUyCu[1]+(27435.68479189101*rdxLy[1]+9976.61265159673*rdxUx[0]+9976.61265159673*rdxLx[0])*rdxUySq[1]+((-27435.68479189101*rdxLySq[1])+9560.920457780201*rdxUxSq[0]-7898.15168251408*rdxLx[0]*rdxUx[0]+9560.920457780201*rdxLxSq[0])*rdxUy[1]-415.6921938165305*rdxLyCu[1]+((-9976.61265159673*rdxUx[0])-9976.61265159673*rdxLx[0])*rdxLySq[1]+((-9560.920457780201*rdxUxSq[0])+7898.15168251408*rdxLx[0]*rdxUx[0]-9560.920457780201*rdxLxSq[0])*rdxLy[1])*rho[3]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[2]+((-1104.0*rdxUyCu[1])+((-75072.0*rdxLy[1])-27600.0*rdxUx[0]-27600.0*rdxLx[0])*rdxUySq[1]+((-75072.0*rdxLySq[1])+((-126960.0*rdxUx[0])-126960.0*rdxLx[0])*rdxLy[1]-26928.0*rdxUxSq[0]-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxUy[1]-1104.0*rdxLyCu[1]+((-27600.0*rdxUx[0])-27600.0*rdxLx[0])*rdxLySq[1]+((-26928.0*rdxUxSq[0])-81936.0*rdxLx[0]*rdxUx[0]-26928.0*rdxLxSq[0])*rdxLy[1]-432.0*rdxUxCu[0]-29376.0*rdxLx[0]*rdxUxSq[0]-29376.0*rdxLxSq[0]*rdxUx[0]-432.0*rdxLxCu[0])*rho[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxUySq[1]+((60164.51685171252*rdxUx[0]-60164.51685171252*rdxLx[0])*rho[0]*rdxLy[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0])*rdxUy[1]+(65208.24880335309*rdxUx[0]-65208.24880335309*rdxLx[0])*rho[0]*rdxLySq[1]+(66289.04850727606*rdxUxSq[0]-66289.04850727606*rdxLxSq[0])*rho[0]*rdxLy[1]+(1080.799703922979*rdxUxCu[0]+71332.78045891662*rdxLx[0]*rdxUxSq[0]-71332.78045891662*rdxLxSq[0]*rdxUx[0]-1080.799703922979*rdxLxCu[0])*rho[0])*omega*volFac+((415.6921938165305*rdxUyR4[1]+(28630.79984911354*rdxLy[1]+10574.170180208*rdxUx[0]+10574.170180208*rdxLx[0])*rdxUyCu[1]+(52637.02404201817*rdxLySq[1]+(68719.1157902952*rdxUx[0]+68719.1157902952*rdxLx[0])*rdxLy[1]+10392.30484541326*rdxUxSq[0]+47804.60228890101*rdxLx[0]*rdxUx[0]+10392.30484541326*rdxLxSq[0])*rdxUySq[1]+(779.4228634059946*rdxLyCu[1]+(19303.70625035514*rdxUx[0]+19303.70625035514*rdxLx[0])*rdxLySq[1]+(18758.11024597094*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0]+18758.11024597094*rdxLxSq[0])*rdxLy[1]+233.8268590217983*rdxUxCu[0]+15900.22641348229*rdxLx[0]*rdxUxSq[0]+15900.22641348229*rdxLxSq[0]*rdxUx[0]+233.8268590217983*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-181.8653347947321*rdxUx[0]*rdxUyCu[1])+((-12003.11209645232*rdxUx[0]*rdxLy[1])+9145.22826396367*rdxUxSq[0]-17874.76433411081*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(12003.11209645232*rdxUx[0]*rdxLySq[1]+9327.093598758403*rdxUxCu[0]+3455.441361099909*rdxLx[0]*rdxUxSq[0]-17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+181.8653347947321*rdxUx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxUxSq[0])*rdxLySq[1]+((-9327.093598758403*rdxUxCu[0])-3455.441361099909*rdxLx[0]*rdxUxSq[0]+17692.89899931608*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((-779.4228634059946*rdxLy[1]*rdxUyCu[1])+(((-19303.70625035514*rdxUx[0])-19303.70625035514*rdxLx[0])*rdxLy[1]-52637.02404201817*rdxLySq[1])*rdxUySq[1]+((-28630.79984911354*rdxLyCu[1])+((-68719.1157902952*rdxUx[0])-68719.1157902952*rdxLx[0])*rdxLySq[1]+((-18758.11024597094*rdxUxSq[0])-40893.71956670118*rdxLx[0]*rdxUx[0]-18758.11024597094*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-415.6921938165305*rdxLyR4[1]+((-10574.170180208*rdxUx[0])-10574.170180208*rdxLx[0])*rdxLyCu[1]+((-10392.30484541326*rdxUxSq[0])-47804.60228890101*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-233.8268590217983*rdxUxCu[0])-15900.22641348229*rdxLx[0]*rdxUxSq[0]-15900.22641348229*rdxLxSq[0]*rdxUx[0]-233.8268590217983*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-181.8653347947321*rdxLx[0]*rdxUyCu[1])+((-12003.11209645232*rdxLx[0]*rdxLy[1])-17874.76433411081*rdxLx[0]*rdxUx[0]+9145.22826396367*rdxLxSq[0])*rdxUySq[1]+(12003.11209645232*rdxLx[0]*rdxLySq[1]-17692.89899931608*rdxLx[0]*rdxUxSq[0]+3455.441361099909*rdxLxSq[0]*rdxUx[0]+9327.093598758403*rdxLxCu[0])*rdxUy[1]+181.8653347947321*rdxLx[0]*rdxLyCu[1]+(17874.76433411081*rdxLx[0]*rdxUx[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+(17692.89899931608*rdxLx[0]*rdxUxSq[0]-3455.441361099909*rdxLxSq[0]*rdxUx[0]-9327.093598758403*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((24375.0*rdxLx[0]-24375.0*rdxUx[0])*rdxUyCu[1]+((32565.0*rdxLx[0]-32565.0*rdxUx[0])*rdxLy[1]-24960.0*rdxUxSq[0]+24960.0*rdxLxSq[0])*rdxUySq[1]+((46215.0*rdxLx[0]-46215.0*rdxUx[0])*rdxLySq[1]+(46800.0*rdxLxSq[0]-46800.0*rdxUxSq[0])*rdxLy[1]-585.0*rdxUxCu[0]-38610.0*rdxLx[0]*rdxUxSq[0]+38610.0*rdxLxSq[0]*rdxUx[0]+585.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(225.0*rdxUx[0]*rdxUyCu[1]+(14850.0*rdxUx[0]*rdxLy[1]-8640.0*rdxUxSq[0]+19440.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-14850.0*rdxUx[0]*rdxLySq[1])-8865.0*rdxUxCu[0]-4275.0*rdxLx[0]*rdxUxSq[0]+19215.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-225.0*rdxUx[0]*rdxLyCu[1]+(8640.0*rdxUxSq[0]-19440.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(8865.0*rdxUxCu[0]+4275.0*rdxLx[0]*rdxUxSq[0]-19215.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+((46215.0*rdxUx[0]-46215.0*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((32565.0*rdxUx[0]-32565.0*rdxLx[0])*rdxLySq[1]+(46800.0*rdxUxSq[0]-46800.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(24375.0*rdxUx[0]-24375.0*rdxLx[0])*rdxLyCu[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLySq[1]+(585.0*rdxUxCu[0]+38610.0*rdxLx[0]*rdxUxSq[0]-38610.0*rdxLxSq[0]*rdxUx[0]-585.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-225.0*rdxLx[0]*rdxUyCu[1])+((-14850.0*rdxLx[0]*rdxLy[1])-19440.0*rdxLx[0]*rdxUx[0]+8640.0*rdxLxSq[0])*rdxUySq[1]+(14850.0*rdxLx[0]*rdxLySq[1]-19215.0*rdxLx[0]*rdxUxSq[0]+4275.0*rdxLxSq[0]*rdxUx[0]+8865.0*rdxLxCu[0])*rdxUy[1]+225.0*rdxLx[0]*rdxLyCu[1]+(19440.0*rdxLx[0]*rdxUx[0]-8640.0*rdxLxSq[0])*rdxLySq[1]+(19215.0*rdxLx[0]*rdxUxSq[0]-4275.0*rdxLxSq[0]*rdxUx[0]-8865.0*rdxLxCu[0])*rdxLy[1])*phiLx[2]+(36.0*phiC[1]-396.0*phiUy[1])*rdxUyR4[1]+(((-27378.0*phiUy[1])-846.0*phiLy[1]+4824.0*phiC[1])*rdxLy[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiUy[1]+483.0*rdxUx[0]*phiUx[1]+483.0*rdxLx[0]*phiLx[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*phiC[1]+(23169.64365284887*phiUy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiUy[0])*rdxLx[0])*rdxUyCu[1]+(((-57078.0*phiUy[1])-57078.0*phiLy[1]+161676.0*phiC[1])*rdxLySq[1]+(((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*phiC[1]+(33842.54072908828*phiUy[0]-40633.91194556586*phiUx[0]+50189.63625092335*phiLy[0])*rdxUx[0]+((-33842.54072908828*phiUy[0])-50189.63625092335*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLy[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiUy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*phiC[1]+(23777.59348630554*phiUy[0]+21740.70173660454*phiUx[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiUy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-846.0*phiUy[1])-27378.0*phiLy[1]+4824.0*phiC[1])*rdxLyCu[1]+(((-20925.0*rdxUx[0])-20925.0*rdxLx[0])*phiUy[1]+32844.0*rdxUx[0]*phiUx[1]+((-71415.0*rdxUx[0])-71415.0*rdxLx[0])*phiLy[1]+32844.0*rdxLx[0]*phiLx[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*phiC[1]+(50189.63625092335*phiUy[0]-40633.91194556586*phiUx[0]+33842.54072908828*phiLy[0])*rdxUx[0]+((-50189.63625092335*phiUy[0])-33842.54072908828*phiLy[0]+40633.91194556586*phiLx[0])*rdxLx[0])*rdxLySq[1]+(((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiUy[1]+(22980.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0])*phiUx[1]+((-20322.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0]-20322.0*rdxLxSq[0])*phiLy[1]+(88110.0*rdxLx[0]*rdxUx[0]+22980.0*rdxLxSq[0])*phiLx[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*phiC[1]+(50797.58608438003*phiUy[0]-34876.57506120691*phiUx[0]+50797.58608438003*phiLy[0])*rdxUxSq[0]+(102561.6565193835*phiLx[0]-102561.6565193835*phiUx[0])*rdxLx[0]*rdxUx[0]+((-50797.58608438003*phiUy[0])-50797.58608438003*phiLy[0]+34876.57506120691*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiUy[1]+((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*phiC[1]+(607.9498334566757*phiUy[0]+22712.38223965068*phiUx[0])*rdxUxCu[0]+(40124.68900814059*phiUy[0]-44349.16092780109*phiUx[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40124.68900814059*phiUy[0])-51862.79733103487*phiUx[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiUy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxUy[1]+(36.0*phiC[1]-396.0*phiLy[1])*rdxLyR4[1]+(483.0*rdxUx[0]*phiUx[1]+((-10125.0*rdxUx[0])-10125.0*rdxLx[0])*phiLy[1]+483.0*rdxLx[0]*phiLx[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*phiC[1]+(23169.64365284887*phiLy[0]-597.5575286112626*phiUx[0])*rdxUx[0]+(597.5575286112626*phiLx[0]-23169.64365284887*phiLy[0])*rdxLx[0])*rdxLyCu[1]+((47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxUxSq[0])*phiUx[1]+((-9972.0*rdxUxSq[0])-50364.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*phiLy[1]+(47370.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*phiLx[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*phiC[1]+(21740.70173660454*phiUx[0]+23777.59348630554*phiLy[0])*rdxUxSq[0]+(51618.57816716767*phiLx[0]-51618.57816716767*phiUx[0])*rdxLx[0]*rdxUx[0]+((-23777.59348630554*phiLy[0])-21740.70173660454*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-24099.0*rdxUxCu[0])+35847.0*rdxLx[0]*rdxUxSq[0]+47661.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-243.0*rdxUxCu[0])-16524.0*rdxLx[0]*rdxUxSq[0]-16524.0*rdxLxSq[0]*rdxUx[0]-243.0*rdxLxCu[0])*phiLy[1]+(47661.0*rdxLx[0]*rdxUxSq[0]+35847.0*rdxLxSq[0]*rdxUx[0]-24099.0*rdxLxCu[0])*phiLx[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*phiC[1]+(22712.38223965068*phiUx[0]+607.9498334566757*phiLy[0])*rdxUxCu[0]+((-44349.16092780109*phiUx[0])+40124.68900814059*phiLy[0]+51862.79733103487*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-51862.79733103487*phiUx[0])-40124.68900814059*phiLy[0]+44349.16092780109*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-607.9498334566757*phiLy[0])-22712.38223965068*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-396.0*rdxUxR4[0])-25758.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]+774.0*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(774.0*rdxLx[0]*rdxUxCu[0]+51462.0*rdxLxSq[0]*rdxUxSq[0]-25758.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[1]+(36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0])*phiC[1]+374.1229744348773*phiUx[0]*rdxUxR4[0]+(24224.46259465831*phiUx[0]+841.7766924784738*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(56024.91542162288*phiLx[0]-56024.91542162288*phiUx[0])*rdxLxSq[0]*rdxUxSq[0]+((-841.7766924784738*phiUx[0])-24224.46259465831*phiLx[0])*rdxLxCu[0]*rdxUx[0]-374.1229744348773*phiLx[0]*rdxLxR4[0])*omega-36.0*phiC[1]*rdxUyR4[1]+(((-3072.0*rdxUx[0])-3072.0*rdxLx[0])*phiC[1]-4824.0*phiC[1]*rdxLy[1])*rdxUyCu[1]+((-161676.0*phiC[1]*rdxLySq[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*phiC[1]*rdxLy[1]+((-6072.0*rdxUxSq[0])-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*phiC[1])*rdxUySq[1]+((-4824.0*phiC[1]*rdxLyCu[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*phiC[1]*rdxLySq[1]+((-207144.0*rdxUxSq[0])-500088.0*rdxLx[0]*rdxUx[0]-207144.0*rdxLxSq[0])*phiC[1]*rdxLy[1]+((-3072.0*rdxUxCu[0])-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*phiC[1])*rdxUy[1]-36.0*phiC[1]*rdxLyR4[1]+((-3072.0*rdxUx[0])-3072.0*rdxLx[0])*phiC[1]*rdxLyCu[1]+((-6072.0*rdxUxSq[0])-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*phiC[1]*rdxLySq[1]+((-3072.0*rdxUxCu[0])-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*phiC[1]*rdxLy[1]+((-36.0*rdxUxR4[0])-4824.0*rdxLx[0]*rdxUxCu[0]-161676.0*rdxLxSq[0]*rdxUxSq[0]-4824.0*rdxLxCu[0]*rdxUx[0]-36.0*rdxLxR4[0])*phiC[1]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[2] = -(1.0*((((9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxUySq[1]+((7898.15168251408*rdxLx[0]-7898.15168251408*rdxUx[0])*rdxLy[1]+9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxUy[1]+(9560.920457780201*rdxUx[0]-9560.920457780201*rdxLx[0])*rdxLySq[1]+(9976.61265159673*rdxUxSq[0]-9976.61265159673*rdxLxSq[0])*rdxLy[1]+415.6921938165305*rdxUxCu[0]+27435.68479189101*rdxLx[0]*rdxUxSq[0]-27435.68479189101*rdxLxSq[0]*rdxUx[0]-415.6921938165305*rdxLxCu[0])*rho[3]+((-432.0*rdxUyCu[1])+((-29376.0*rdxLy[1])-26928.0*rdxUx[0]-26928.0*rdxLx[0])*rdxUySq[1]+((-29376.0*rdxLySq[1])+((-81936.0*rdxUx[0])-81936.0*rdxLx[0])*rdxLy[1]-27600.0*rdxUxSq[0]-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxUy[1]-432.0*rdxLyCu[1]+((-26928.0*rdxUx[0])-26928.0*rdxLx[0])*rdxLySq[1]+((-27600.0*rdxUxSq[0])-126960.0*rdxLx[0]*rdxUx[0]-27600.0*rdxLxSq[0])*rdxLy[1]-1104.0*rdxUxCu[0]-75072.0*rdxLx[0]*rdxUxSq[0]-75072.0*rdxLxSq[0]*rdxUx[0]-1104.0*rdxLxCu[0])*rho[2]+((24960.0*rdxLx[0]-24960.0*rdxUx[0])*rdxUySq[1]+(24960.0*rdxLxSq[0]-24960.0*rdxUxSq[0])*rdxUy[1]+(24960.0*rdxUx[0]-24960.0*rdxLx[0])*rdxLySq[1]+(24960.0*rdxUxSq[0]-24960.0*rdxLxSq[0])*rdxLy[1])*rho[1]+1080.799703922979*rho[0]*rdxUyCu[1]+(71332.78045891662*rho[0]*rdxLy[1]+(66289.04850727606*rdxUx[0]+66289.04850727606*rdxLx[0])*rho[0])*rdxUySq[1]+((65208.24880335309*rdxUxSq[0]+60164.51685171252*rdxLx[0]*rdxUx[0]+65208.24880335309*rdxLxSq[0])*rho[0]-71332.78045891662*rho[0]*rdxLySq[1])*rdxUy[1]-1080.799703922979*rho[0]*rdxLyCu[1]+((-66289.04850727606*rdxUx[0])-66289.04850727606*rdxLx[0])*rho[0]*rdxLySq[1]+((-65208.24880335309*rdxUxSq[0])-60164.51685171252*rdxLx[0]*rdxUx[0]-65208.24880335309*rdxLxSq[0])*rho[0]*rdxLy[1])*omega*volFac+(((9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxUyCu[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLy[1]+9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxUySq[1]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1]-181.8653347947321*rdxUxCu[0]-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(233.8268590217983*rdxUx[0]*rdxUyCu[1]+(15900.22641348229*rdxUx[0]*rdxLy[1]+10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(15900.22641348229*rdxUx[0]*rdxLySq[1]+(47804.60228890101*rdxUxSq[0]+40893.71956670118*rdxLx[0]*rdxUx[0])*rdxLy[1]+10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+233.8268590217983*rdxUx[0]*rdxLyCu[1]+(10392.30484541326*rdxUxSq[0]+18758.11024597094*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(10574.170180208*rdxUxCu[0]+68719.1157902952*rdxLx[0]*rdxUxSq[0]+19303.70625035514*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+415.6921938165305*rdxUxR4[0]+28630.79984911354*rdxLx[0]*rdxUxCu[0]+52637.02404201817*rdxLxSq[0]*rdxUxSq[0]+779.4228634059946*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((17692.89899931608*rdxLx[0]-17692.89899931608*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((3455.441361099909*rdxUx[0]-3455.441361099909*rdxLx[0])*rdxLySq[1]+(17874.76433411081*rdxLxSq[0]-17874.76433411081*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(9327.093598758403*rdxUx[0]-9327.093598758403*rdxLx[0])*rdxLyCu[1]+(9145.22826396367*rdxUxSq[0]-9145.22826396367*rdxLxSq[0])*rdxLySq[1]+((-181.8653347947321*rdxUxCu[0])-12003.11209645232*rdxLx[0]*rdxUxSq[0]+12003.11209645232*rdxLxSq[0]*rdxUx[0]+181.8653347947321*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-233.8268590217983*rdxLx[0]*rdxUyCu[1])+((-15900.22641348229*rdxLx[0]*rdxLy[1])-18758.11024597094*rdxLx[0]*rdxUx[0]-10392.30484541326*rdxLxSq[0])*rdxUySq[1]+((-15900.22641348229*rdxLx[0]*rdxLySq[1])+((-40893.71956670118*rdxLx[0]*rdxUx[0])-47804.60228890101*rdxLxSq[0])*rdxLy[1]-19303.70625035514*rdxLx[0]*rdxUxSq[0]-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxUy[1]-233.8268590217983*rdxLx[0]*rdxLyCu[1]+((-18758.11024597094*rdxLx[0]*rdxUx[0])-10392.30484541326*rdxLxSq[0])*rdxLySq[1]+((-19303.70625035514*rdxLx[0]*rdxUxSq[0])-68719.1157902952*rdxLxSq[0]*rdxUx[0]-10574.170180208*rdxLxCu[0])*rdxLy[1]-779.4228634059946*rdxLx[0]*rdxUxCu[0]-52637.02404201817*rdxLxSq[0]*rdxUxSq[0]-28630.79984911354*rdxLxCu[0]*rdxUx[0]-415.6921938165305*rdxLxR4[0])*phiLx[3]+((-396.0*rdxUyR4[1])+((-25758.0*rdxLy[1])-24099.0*rdxUx[0]-24099.0*rdxLx[0])*rdxUyCu[1]+(51462.0*rdxLySq[1]+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLy[1]-23220.0*rdxUxSq[0]+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxUySq[1]+(774.0*rdxLyCu[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1]+483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-243.0*rdxUx[0]*rdxUyCu[1])+((-16524.0*rdxUx[0]*rdxLy[1])-9972.0*rdxUxSq[0]-20322.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-16524.0*rdxUx[0]*rdxLySq[1])+((-50364.0*rdxUxSq[0])-41814.0*rdxLx[0]*rdxUx[0])*rdxLy[1]-10125.0*rdxUxCu[0]-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-243.0*rdxUx[0]*rdxLyCu[1]+((-9972.0*rdxUxSq[0])-20322.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-10125.0*rdxUxCu[0])-71415.0*rdxLx[0]*rdxUxSq[0]-20925.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-396.0*rdxUxR4[0]-27378.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-846.0*rdxLxCu[0]*rdxUx[0])*phiUx[2]+(774.0*rdxLy[1]*rdxUyCu[1]+(51462.0*rdxLySq[1]+(47661.0*rdxUx[0]+47661.0*rdxLx[0])*rdxLy[1])*rdxUySq[1]+((-25758.0*rdxLyCu[1])+(35847.0*rdxUx[0]+35847.0*rdxLx[0])*rdxLySq[1]+(47370.0*rdxUxSq[0]+88110.0*rdxLx[0]*rdxUx[0]+47370.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-396.0*rdxLyR4[1]+((-24099.0*rdxUx[0])-24099.0*rdxLx[0])*rdxLyCu[1]+((-23220.0*rdxUxSq[0])+22980.0*rdxLx[0]*rdxUx[0]-23220.0*rdxLxSq[0])*rdxLySq[1]+(483.0*rdxUxCu[0]+32844.0*rdxLx[0]*rdxUxSq[0]+32844.0*rdxLxSq[0]*rdxUx[0]+483.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-243.0*rdxLx[0]*rdxUyCu[1])+((-16524.0*rdxLx[0]*rdxLy[1])-20322.0*rdxLx[0]*rdxUx[0]-9972.0*rdxLxSq[0])*rdxUySq[1]+((-16524.0*rdxLx[0]*rdxLySq[1])+((-41814.0*rdxLx[0]*rdxUx[0])-50364.0*rdxLxSq[0])*rdxLy[1]-20925.0*rdxLx[0]*rdxUxSq[0]-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxUy[1]-243.0*rdxLx[0]*rdxLyCu[1]+((-20322.0*rdxLx[0]*rdxUx[0])-9972.0*rdxLxSq[0])*rdxLySq[1]+((-20925.0*rdxLx[0]*rdxUxSq[0])-71415.0*rdxLxSq[0]*rdxUx[0]-10125.0*rdxLxCu[0])*rdxLy[1]-846.0*rdxLx[0]*rdxUxCu[0]-57078.0*rdxLxSq[0]*rdxUxSq[0]-27378.0*rdxLxCu[0]*rdxUx[0]-396.0*rdxLxR4[0])*phiLx[2]+(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0])*phiC[2]+374.1229744348773*phiUy[0]*rdxUyR4[1]+((24224.46259465831*phiUy[0]+841.7766924784738*phiLy[0])*rdxLy[1]+(8865.0*rdxLx[0]-8865.0*rdxUx[0])*phiUy[1]-585.0*rdxUx[0]*phiUx[1]+585.0*rdxLx[0]*phiLx[1]+(22712.38223965068*phiUy[0]+607.9498334566757*phiUx[0])*rdxUx[0]+(22712.38223965068*phiUy[0]+607.9498334566757*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((56024.91542162288*phiLy[0]-56024.91542162288*phiUy[0])*rdxLySq[1]+((4275.0*rdxLx[0]-4275.0*rdxUx[0])*phiUy[1]-38610.0*rdxUx[0]*phiUx[1]+(19215.0*rdxLx[0]-19215.0*rdxUx[0])*phiLy[1]+38610.0*rdxLx[0]*phiLx[1]+((-44349.16092780109*phiUy[0])+40124.68900814059*phiUx[0]+51862.79733103487*phiLy[0])*rdxUx[0]+((-44349.16092780109*phiUy[0])+51862.79733103487*phiLy[0]+40124.68900814059*phiLx[0])*rdxLx[0])*rdxLy[1]+(8640.0*rdxLxSq[0]-8640.0*rdxUxSq[0])*phiUy[1]+((-24960.0*rdxUxSq[0])-46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(46800.0*rdxLx[0]*rdxUx[0]+24960.0*rdxLxSq[0])*phiLx[1]+(21740.70173660454*phiUy[0]+23777.59348630554*phiUx[0])*rdxUxSq[0]+((-34876.57506120691*phiUy[0])+50797.58608438003*phiUx[0]+50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+(21740.70173660454*phiUy[0]+23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-841.7766924784738*phiUy[0])-24224.46259465831*phiLy[0])*rdxLyCu[1]+((19215.0*rdxUx[0]-19215.0*rdxLx[0])*phiUy[1]+38610.0*rdxUx[0]*phiUx[1]+(4275.0*rdxUx[0]-4275.0*rdxLx[0])*phiLy[1]-38610.0*rdxLx[0]*phiLx[1]+((-51862.79733103487*phiUy[0])-40124.68900814059*phiUx[0]+44349.16092780109*phiLy[0])*rdxUx[0]+((-51862.79733103487*phiUy[0])+44349.16092780109*phiLy[0]-40124.68900814059*phiLx[0])*rdxLx[0])*rdxLySq[1]+((19440.0*rdxUxSq[0]-19440.0*rdxLxSq[0])*phiUy[1]+(19440.0*rdxLxSq[0]-19440.0*rdxUxSq[0])*phiLy[1]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxUxSq[0]+(102561.6565193835*phiLy[0]-102561.6565193835*phiUy[0])*rdxLx[0]*rdxUx[0]+(51618.57816716767*phiLy[0]-51618.57816716767*phiUy[0])*rdxLxSq[0])*rdxLy[1]+(225.0*rdxUxCu[0]+14850.0*rdxLx[0]*rdxUxSq[0]-14850.0*rdxLxSq[0]*rdxUx[0]-225.0*rdxLxCu[0])*phiUy[1]+((-24375.0*rdxUxCu[0])-32565.0*rdxLx[0]*rdxUxSq[0]-46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(46215.0*rdxLx[0]*rdxUxSq[0]+32565.0*rdxLxSq[0]*rdxUx[0]+24375.0*rdxLxCu[0])*phiLx[1]+(23169.64365284887*phiUx[0]-597.5575286112626*phiUy[0])*rdxUxCu[0]+((-40633.91194556586*phiUy[0])+33842.54072908828*phiUx[0]+50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-40633.91194556586*phiUy[0])+50189.63625092335*phiUx[0]+33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(23169.64365284887*phiLx[0]-597.5575286112626*phiUy[0])*rdxLxCu[0])*rdxUy[1]-374.1229744348773*phiLy[0]*rdxLyR4[1]+(585.0*rdxUx[0]*phiUx[1]+(8865.0*rdxUx[0]-8865.0*rdxLx[0])*phiLy[1]-585.0*rdxLx[0]*phiLx[1]+((-607.9498334566757*phiUx[0])-22712.38223965068*phiLy[0])*rdxUx[0]+((-22712.38223965068*phiLy[0])-607.9498334566757*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((24960.0*rdxUxSq[0]+46800.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(8640.0*rdxUxSq[0]-8640.0*rdxLxSq[0])*phiLy[1]+((-46800.0*rdxLx[0]*rdxUx[0])-24960.0*rdxLxSq[0])*phiLx[1]+((-23777.59348630554*phiUx[0])-21740.70173660454*phiLy[0])*rdxUxSq[0]+((-50797.58608438003*phiUx[0])+34876.57506120691*phiLy[0]-50797.58608438003*phiLx[0])*rdxLx[0]*rdxUx[0]+((-21740.70173660454*phiLy[0])-23777.59348630554*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((24375.0*rdxUxCu[0]+32565.0*rdxLx[0]*rdxUxSq[0]+46215.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-225.0*rdxUxCu[0])-14850.0*rdxLx[0]*rdxUxSq[0]+14850.0*rdxLxSq[0]*rdxUx[0]+225.0*rdxLxCu[0])*phiLy[1]+((-46215.0*rdxLx[0]*rdxUxSq[0])-32565.0*rdxLxSq[0]*rdxUx[0]-24375.0*rdxLxCu[0])*phiLx[1]+(597.5575286112626*phiLy[0]-23169.64365284887*phiUx[0])*rdxUxCu[0]+((-33842.54072908828*phiUx[0])+40633.91194556586*phiLy[0]-50189.63625092335*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-50189.63625092335*phiUx[0])+40633.91194556586*phiLy[0]-33842.54072908828*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(597.5575286112626*phiLy[0]-23169.64365284887*phiLx[0])*rdxLxCu[0])*rdxLy[1])*omega+((-36.0*rdxUyR4[1])+((-4824.0*rdxLy[1])-3072.0*rdxUx[0]-3072.0*rdxLx[0])*rdxUyCu[1]+((-161676.0*rdxLySq[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*rdxLy[1]-6072.0*rdxUxSq[0]-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*rdxUySq[1]+((-4824.0*rdxLyCu[1])+((-208896.0*rdxUx[0])-208896.0*rdxLx[0])*rdxLySq[1]+((-207144.0*rdxUxSq[0])-500088.0*rdxLx[0]*rdxUx[0]-207144.0*rdxLxSq[0])*rdxLy[1]-3072.0*rdxUxCu[0]-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*rdxUy[1]-36.0*rdxLyR4[1]+((-3072.0*rdxUx[0])-3072.0*rdxLx[0])*rdxLyCu[1]+((-6072.0*rdxUxSq[0])-207144.0*rdxLx[0]*rdxUx[0]-6072.0*rdxLxSq[0])*rdxLySq[1]+((-3072.0*rdxUxCu[0])-208896.0*rdxLx[0]*rdxUxSq[0]-208896.0*rdxLxSq[0]*rdxUx[0]-3072.0*rdxLxCu[0])*rdxLy[1]-36.0*rdxUxR4[0]-4824.0*rdxLx[0]*rdxUxCu[0]-161676.0*rdxLxSq[0]*rdxUxSq[0]-4824.0*rdxLxCu[0]*rdxUx[0]-36.0*rdxLxR4[0])*phiC[2]))/(36.0*rdxUyR4[1]+(4824.0*rdxLy[1]+3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxUyCu[1]+(161676.0*rdxLySq[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLy[1]+6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxUySq[1]+(4824.0*rdxLyCu[1]+(208896.0*rdxUx[0]+208896.0*rdxLx[0])*rdxLySq[1]+(207144.0*rdxUxSq[0]+500088.0*rdxLx[0]*rdxUx[0]+207144.0*rdxLxSq[0])*rdxLy[1]+3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxUy[1]+36.0*rdxLyR4[1]+(3072.0*rdxUx[0]+3072.0*rdxLx[0])*rdxLyCu[1]+(6072.0*rdxUxSq[0]+207144.0*rdxLx[0]*rdxUx[0]+6072.0*rdxLxSq[0])*rdxLySq[1]+(3072.0*rdxUxCu[0]+208896.0*rdxLx[0]*rdxUxSq[0]+208896.0*rdxLxSq[0]*rdxUx[0]+3072.0*rdxLxCu[0])*rdxLy[1]+36.0*rdxUxR4[0]+4824.0*rdxLx[0]*rdxUxCu[0]+161676.0*rdxLxSq[0]*rdxUxSq[0]+4824.0*rdxLxCu[0]*rdxUx[0]+36.0*rdxLxR4[0]); 
  phiC[3] = (((144.0*rdxUyCu[1]+(9792.0*rdxLy[1]+3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxUySq[1]+(9792.0*rdxLySq[1]+(31568.0*rdxUx[0]+31568.0*rdxLx[0])*rdxLy[1]+3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxUy[1]+144.0*rdxLyCu[1]+(3824.0*rdxUx[0]+3824.0*rdxLx[0])*rdxLySq[1]+(3824.0*rdxUxSq[0]+31568.0*rdxLx[0]*rdxUx[0]+3824.0*rdxLxSq[0])*rdxLy[1]+144.0*rdxUxCu[0]+9792.0*rdxLx[0]*rdxUxSq[0]+9792.0*rdxLxSq[0]*rdxUx[0]+144.0*rdxLxCu[0])*rho[3]+((8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxUySq[1]+((6845.064791512203*rdxUx[0]-6845.064791512203*rdxLx[0])*rdxLy[1]-8646.397631383834*rdxUxSq[0]+8646.397631383834*rdxLxSq[0])*rdxUy[1]+(8286.131063409508*rdxLx[0]-8286.131063409508*rdxUx[0])*rdxLySq[1]+(8646.397631383834*rdxLxSq[0]-8646.397631383834*rdxUxSq[0])*rdxLy[1]-360.2665679743264*rdxUxCu[0]-23777.59348630554*rdxLx[0]*rdxUxSq[0]+23777.59348630554*rdxLxSq[0]*rdxUx[0]+360.2665679743264*rdxLxCu[0])*rho[2]+((-360.2665679743264*rdxUyCu[1])+((-23777.59348630554*rdxLy[1])-8646.397631383834*rdxUx[0]-8646.397631383834*rdxLx[0])*rdxUySq[1]+(23777.59348630554*rdxLySq[1]-8286.131063409508*rdxUxSq[0]+6845.064791512203*rdxLx[0]*rdxUx[0]-8286.131063409508*rdxLxSq[0])*rdxUy[1]+360.2665679743264*rdxLyCu[1]+(8646.397631383834*rdxUx[0]+8646.397631383834*rdxLx[0])*rdxLySq[1]+(8286.131063409508*rdxUxSq[0]-6845.064791512203*rdxLx[0]*rdxUx[0]+8286.131063409508*rdxLxSq[0])*rdxLy[1])*rho[1]+(21632.0*rdxUx[0]-21632.0*rdxLx[0])*rho[0]*rdxUySq[1]+(21632.0*rdxUxSq[0]-21632.0*rdxLxSq[0])*rho[0]*rdxUy[1]+(21632.0*rdxLx[0]-21632.0*rdxUx[0])*rho[0]*rdxLySq[1]+(21632.0*rdxLxSq[0]-21632.0*rdxUxSq[0])*rho[0]*rdxLy[1])*omega*volFac+((132.0*rdxUyR4[1]+(8586.0*rdxLy[1]+3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxUyCu[1]+((-17154.0*rdxLySq[1])+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLy[1]+2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-258.0*rdxLyCu[1])+((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1]-63.0*rdxUxCu[0]-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-63.0*rdxUx[0]*rdxUyCu[1])+((-4284.0*rdxUx[0]*rdxLy[1])+2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-4284.0*rdxUx[0]*rdxLySq[1])+((-17516.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0])*rdxLy[1]+3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-63.0*rdxUx[0]*rdxLyCu[1]+(2812.0*rdxUxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(3007.0*rdxUxCu[0]-13811.0*rdxLx[0]*rdxUxSq[0]-6353.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+132.0*rdxUxR4[0]+8586.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]-258.0*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((-258.0*rdxLy[1]*rdxUyCu[1])+(((-6353.0*rdxUx[0])-6353.0*rdxLx[0])*rdxLy[1]-17154.0*rdxLySq[1])*rdxUySq[1]+(8586.0*rdxLyCu[1]+((-13811.0*rdxUx[0])-13811.0*rdxLx[0])*rdxLySq[1]+((-6158.0*rdxUxSq[0])-10106.0*rdxLx[0]*rdxUx[0]-6158.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+132.0*rdxLyR4[1]+(3007.0*rdxUx[0]+3007.0*rdxLx[0])*rdxLyCu[1]+(2812.0*rdxUxSq[0]-17516.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxLySq[1]+((-63.0*rdxUxCu[0])-4284.0*rdxLx[0]*rdxUxSq[0]-4284.0*rdxLxSq[0]*rdxUx[0]-63.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-63.0*rdxLx[0]*rdxUyCu[1])+((-4284.0*rdxLx[0]*rdxLy[1])-6158.0*rdxLx[0]*rdxUx[0]+2812.0*rdxLxSq[0])*rdxUySq[1]+((-4284.0*rdxLx[0]*rdxLySq[1])+((-10106.0*rdxLx[0]*rdxUx[0])-17516.0*rdxLxSq[0])*rdxLy[1]-6353.0*rdxLx[0]*rdxUxSq[0]-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxUy[1]-63.0*rdxLx[0]*rdxLyCu[1]+(2812.0*rdxLxSq[0]-6158.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6353.0*rdxLx[0]*rdxUxSq[0])-13811.0*rdxLxSq[0]*rdxUx[0]+3007.0*rdxLxCu[0])*rdxLy[1]-258.0*rdxLx[0]*rdxUxCu[0]-17154.0*rdxLxSq[0]*rdxUxSq[0]+8586.0*rdxLxCu[0]*rdxUx[0]+132.0*rdxLxR4[0])*phiLx[3]+((-12.0*rdxUyR4[1])+((-1608.0*rdxLy[1])-1024.0*rdxUx[0]-1024.0*rdxLx[0])*rdxUyCu[1]+((-53892.0*rdxLySq[1])+((-69632.0*rdxUx[0])-69632.0*rdxLx[0])*rdxLy[1]-2024.0*rdxUxSq[0]-69048.0*rdxLx[0]*rdxUx[0]-2024.0*rdxLxSq[0])*rdxUySq[1]+((-1608.0*rdxLyCu[1])+((-69632.0*rdxUx[0])-69632.0*rdxLx[0])*rdxLySq[1]+((-69048.0*rdxUxSq[0])-166696.0*rdxLx[0]*rdxUx[0]-69048.0*rdxLxSq[0])*rdxLy[1]-1024.0*rdxUxCu[0]-69632.0*rdxLx[0]*rdxUxSq[0]-69632.0*rdxLxSq[0]*rdxUx[0]-1024.0*rdxLxCu[0])*rdxUy[1]-12.0*rdxLyR4[1]+((-1024.0*rdxUx[0])-1024.0*rdxLx[0])*rdxLyCu[1]+((-2024.0*rdxUxSq[0])-69048.0*rdxLx[0]*rdxUx[0]-2024.0*rdxLxSq[0])*rdxLySq[1]+((-1024.0*rdxUxCu[0])-69632.0*rdxLx[0]*rdxUxSq[0]-69632.0*rdxLxSq[0]*rdxUx[0]-1024.0*rdxLxCu[0])*rdxLy[1]-12.0*rdxUxR4[0]-1608.0*rdxLx[0]*rdxUxCu[0]-53892.0*rdxLxSq[0]*rdxUxSq[0]-1608.0*rdxLxCu[0]*rdxUx[0]-12.0*rdxLxR4[0])*phiC[3]+((8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxUyCu[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLy[1]-7925.864495435182*rdxUxSq[0]+7925.864495435182*rdxLxSq[0])*rdxUySq[1]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1]+157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(77.94228634059945*rdxUx[0]*rdxUyCu[1]+(5300.075471160763*rdxUx[0]*rdxLy[1]-2591.14800812304*rdxUxSq[0]+6730.749438212657*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(5300.075471160763*rdxUx[0]*rdxLySq[1]+(20937.03016189259*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0])*rdxLy[1]-2793.797952608599*rdxUxCu[0]+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+77.94228634059945*rdxUx[0]*rdxLyCu[1]+(6730.749438212657*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxUxSq[0])*rdxLySq[1]+((-2793.797952608599*rdxUxCu[0])+17086.68121666697*rdxLx[0]*rdxUxSq[0]+6933.399382698215*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-124.7076581449591*rdxUxR4[0]-8074.820864886104*rdxLx[0]*rdxUxCu[0]+18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+280.592230826158*rdxLxCu[0]*rdxUx[0])*phiUx[2]+((15333.84579940727*rdxUx[0]-15333.84579940727*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((2994.715846286589*rdxLx[0]-2994.715846286589*rdxUx[0])*rdxLySq[1]+(15491.46242289604*rdxUxSq[0]-15491.46242289604*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(8083.48111892395*rdxLx[0]-8083.48111892395*rdxUx[0])*rdxLyCu[1]+(7925.864495435182*rdxLxSq[0]-7925.864495435182*rdxUxSq[0])*rdxLySq[1]+(157.6166234887678*rdxUxCu[0]+10402.69715025867*rdxLx[0]*rdxUxSq[0]-10402.69715025867*rdxLxSq[0]*rdxUx[0]-157.6166234887678*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-77.94228634059945*rdxLx[0]*rdxUyCu[1])+((-5300.075471160763*rdxLx[0]*rdxLy[1])-6730.749438212657*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*rdxUySq[1]+((-5300.075471160763*rdxLx[0]*rdxLySq[1])+((-13236.33227144136*rdxLx[0]*rdxUx[0])-20937.03016189259*rdxLxSq[0])*rdxLy[1]-6933.399382698215*rdxLx[0]*rdxUxSq[0]-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxUy[1]-77.94228634059945*rdxLx[0]*rdxLyCu[1]+(2591.14800812304*rdxLxSq[0]-6730.749438212657*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-6933.399382698215*rdxLx[0]*rdxUxSq[0])-17086.68121666697*rdxLxSq[0]*rdxUx[0]+2793.797952608599*rdxLxCu[0])*rdxLy[1]-280.592230826158*rdxLx[0]*rdxUxCu[0]-18674.97180720763*rdxLxSq[0]*rdxUxSq[0]+8074.820864886104*rdxLxCu[0]*rdxUx[0]+124.7076581449591*rdxLxR4[0])*phiLx[2]-124.7076581449591*phiUy[1]*rdxUyR4[1]+(((-8074.820864886104*phiUy[1])-280.592230826158*phiLy[1])*rdxLy[1]+((-2793.797952608599*rdxUx[0])-2793.797952608599*rdxLx[0])*phiUy[1]+157.6166234887678*rdxUx[0]*phiUx[1]+157.6166234887678*rdxLx[0]*phiLx[1]+(7683.0*phiUy[0]-195.0*phiUx[0])*rdxUx[0]+(195.0*phiLx[0]-7683.0*phiUy[0])*rdxLx[0])*rdxUyCu[1]+((18674.97180720763*phiUy[1]-18674.97180720763*phiLy[1])*rdxLySq[1]+((17086.68121666697*rdxUx[0]+17086.68121666697*rdxLx[0])*phiUy[1]+10402.69715025867*rdxUx[0]*phiUx[1]+((-6933.399382698215*rdxUx[0])-6933.399382698215*rdxLx[0])*phiLy[1]+10402.69715025867*rdxLx[0]*phiLx[1]+(3705.0*phiUy[0]-12870.0*phiUx[0]+16653.0*phiLy[0])*rdxUx[0]+((-3705.0*phiUy[0])-16653.0*phiLy[0]+12870.0*phiLx[0])*rdxLx[0])*rdxLy[1]+((-2591.14800812304*rdxUxSq[0])+20937.03016189259*rdxLx[0]*rdxUx[0]-2591.14800812304*rdxLxSq[0])*phiUy[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxUxSq[0])*phiUx[1]+(15491.46242289604*rdxLx[0]*rdxUx[0]-7925.864495435182*rdxLxSq[0])*phiLx[1]+(7488.0*phiUy[0]+7488.0*phiUx[0])*rdxUxSq[0]+(16848.0*phiLx[0]-16848.0*phiUx[0])*rdxLx[0]*rdxUx[0]+((-7488.0*phiUy[0])-7488.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((280.592230826158*phiUy[1]+8074.820864886104*phiLy[1])*rdxLyCu[1]+((6933.399382698215*rdxUx[0]+6933.399382698215*rdxLx[0])*phiUy[1]-10402.69715025867*rdxUx[0]*phiUx[1]+((-17086.68121666697*rdxUx[0])-17086.68121666697*rdxLx[0])*phiLy[1]-10402.69715025867*rdxLx[0]*phiLx[1]+((-16653.0*phiUy[0])+12870.0*phiUx[0]-3705.0*phiLy[0])*rdxUx[0]+(16653.0*phiUy[0]+3705.0*phiLy[0]-12870.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((6730.749438212657*rdxUxSq[0]+13236.33227144136*rdxLx[0]*rdxUx[0]+6730.749438212657*rdxLxSq[0])*phiUy[1]+((-6730.749438212657*rdxUxSq[0])-13236.33227144136*rdxLx[0]*rdxUx[0]-6730.749438212657*rdxLxSq[0])*phiLy[1]+(16848.0*phiLy[0]-16848.0*phiUy[0])*rdxUxSq[0]+(16848.0*phiUy[0]-16848.0*phiLy[0])*rdxLxSq[0])*rdxLy[1]+(77.94228634059945*rdxUxCu[0]+5300.075471160763*rdxLx[0]*rdxUxSq[0]+5300.075471160763*rdxLxSq[0]*rdxUx[0]+77.94228634059945*rdxLxCu[0])*phiUy[1]+((-8083.48111892395*rdxUxCu[0])-2994.715846286589*rdxLx[0]*rdxUxSq[0]+15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(15333.84579940727*rdxLx[0]*rdxUxSq[0]-2994.715846286589*rdxLxSq[0]*rdxUx[0]-8083.48111892395*rdxLxCu[0])*phiLx[1]+(7683.0*phiUx[0]-195.0*phiUy[0])*rdxUxCu[0]+((-12870.0*phiUy[0])+3705.0*phiUx[0]+16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(12870.0*phiUy[0]-16653.0*phiUx[0]-3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(195.0*phiUy[0]-7683.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+124.7076581449591*phiLy[1]*rdxLyR4[1]+((-157.6166234887678*rdxUx[0]*phiUx[1])+(2793.797952608599*rdxUx[0]+2793.797952608599*rdxLx[0])*phiLy[1]-157.6166234887678*rdxLx[0]*phiLx[1]+(195.0*phiUx[0]-7683.0*phiLy[0])*rdxUx[0]+(7683.0*phiLy[0]-195.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((7925.864495435182*rdxUxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiUx[1]+(2591.14800812304*rdxUxSq[0]-20937.03016189259*rdxLx[0]*rdxUx[0]+2591.14800812304*rdxLxSq[0])*phiLy[1]+(7925.864495435182*rdxLxSq[0]-15491.46242289604*rdxLx[0]*rdxUx[0])*phiLx[1]+((-7488.0*phiUx[0])-7488.0*phiLy[0])*rdxUxSq[0]+(16848.0*phiUx[0]-16848.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(7488.0*phiLy[0]+7488.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((8083.48111892395*rdxUxCu[0]+2994.715846286589*rdxLx[0]*rdxUxSq[0]-15333.84579940727*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-77.94228634059945*rdxUxCu[0])-5300.075471160763*rdxLx[0]*rdxUxSq[0]-5300.075471160763*rdxLxSq[0]*rdxUx[0]-77.94228634059945*rdxLxCu[0])*phiLy[1]+((-15333.84579940727*rdxLx[0]*rdxUxSq[0])+2994.715846286589*rdxLxSq[0]*rdxUx[0]+8083.48111892395*rdxLxCu[0])*phiLx[1]+(195.0*phiLy[0]-7683.0*phiUx[0])*rdxUxCu[0]+((-3705.0*phiUx[0])+12870.0*phiLy[0]-16653.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(16653.0*phiUx[0]-12870.0*phiLy[0]+3705.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(7683.0*phiLx[0]-195.0*phiLy[0])*rdxLxCu[0])*rdxLy[1])*omega+(12.0*rdxUyR4[1]+(1608.0*rdxLy[1]+1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxUyCu[1]+(53892.0*rdxLySq[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLy[1]+2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxUySq[1]+(1608.0*rdxLyCu[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLySq[1]+(69048.0*rdxUxSq[0]+166696.0*rdxLx[0]*rdxUx[0]+69048.0*rdxLxSq[0])*rdxLy[1]+1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxUy[1]+12.0*rdxLyR4[1]+(1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxLyCu[1]+(2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxLySq[1]+(1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxLy[1]+12.0*rdxUxR4[0]+1608.0*rdxLx[0]*rdxUxCu[0]+53892.0*rdxLxSq[0]*rdxUxSq[0]+1608.0*rdxLxCu[0]*rdxUx[0]+12.0*rdxLxR4[0])*phiC[3])/(12.0*rdxUyR4[1]+(1608.0*rdxLy[1]+1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxUyCu[1]+(53892.0*rdxLySq[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLy[1]+2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxUySq[1]+(1608.0*rdxLyCu[1]+(69632.0*rdxUx[0]+69632.0*rdxLx[0])*rdxLySq[1]+(69048.0*rdxUxSq[0]+166696.0*rdxLx[0]*rdxUx[0]+69048.0*rdxLxSq[0])*rdxLy[1]+1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxUy[1]+12.0*rdxLyR4[1]+(1024.0*rdxUx[0]+1024.0*rdxLx[0])*rdxLyCu[1]+(2024.0*rdxUxSq[0]+69048.0*rdxLx[0]*rdxUx[0]+2024.0*rdxLxSq[0])*rdxLySq[1]+(1024.0*rdxUxCu[0]+69632.0*rdxLx[0]*rdxUxSq[0]+69632.0*rdxLxSq[0]*rdxUx[0]+1024.0*rdxLxCu[0])*rdxLy[1]+12.0*rdxUxR4[0]+1608.0*rdxLx[0]*rdxUxCu[0]+53892.0*rdxLxSq[0]*rdxUxSq[0]+1608.0*rdxLxCu[0]*rdxUx[0]+12.0*rdxLxR4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[0]*rho[1]+144.0*rho[0]*rdxCp2[1]+1600.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(81.0*phiUy[0]+81.0*phiLy[0]-162.0*phiC[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+93.53074360871933*rdxCp2[0]*phiUx[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+(900.0*phiUy[0]-54.0*phiUx[0]+900.0*phiLy[0]-2178.0*phiC[0]+864.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiUx[1]+(1020.0*phiUx[0]-2580.0*phiC[0]+3120.0*bcVals[0])*rdxCp2Sq[0])*omega+162.0*phiC[0]*rdxCp2Sq[1]+2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+2580.0*phiC[0]*rdxCp2Sq[0])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]+277.1281292110203*rdxCp2[0]*rho[0])*omega*volFac+(((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1]-162.0*phiC[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]-360.0*rdxCp2[0]*phiUx[1]+189.0*rdxCp2[0]*phiLy[1]-2178.0*rdxCp2[0]*phiC[1]+(155.8845726811989*phiUy[0]+311.7691453623978*phiUx[0]+155.8845726811989*phiLy[0]-1247.076581449591*bcVals[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiUx[1]-2580.0*rdxCp2Sq[0]*phiC[1]+(623.5382907247956*phiUx[0]-1247.076581449591*bcVals[0])*rdxCp2Sq[0])*omega+162.0*phiC[1]*rdxCp2Sq[1]+2178.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+2580.0*rdxCp2Sq[0]*phiC[1])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = ((748.2459488697547*rdxCp2[0]*rho[3]+(368.0*rdxCp2[1]+1600.0*rdxCp2[0])*rho[2])*omega*volFac+((-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3])+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(1020.0*rdxCp2Sq[0]-138.0*rdxCp2[0]*rdxCp2[1])*phiUx[2]+((-161.0*rdxCp2Sq[1])-700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[2]+(199.1858428704209*phiUy[0]-199.1858428704209*phiLy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiUy[0]-866.0254037844386*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[2])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[0]*rho[2])*omega*volFac+(((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiUx[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+623.5382907247956*rdxCp2Sq[0])*phiUx[2]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiUy[0]-150.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[3])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[1]-576.0*rho[0]*rdxCp2[1]-1680.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-324.0*phiUy[0])-324.0*phiLy[0]+648.0*phiC[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+311.7691453623978*rdxCp2[0]*phiUx[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+((-945.0*phiUy[0])-324.0*phiUx[0]-945.0*phiLy[0]+2214.0*phiC[0]+576.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiUx[1]+((-720.0*phiUx[0])+720.0*phiC[0]+2080.0*bcVals[0])*rdxCp2Sq[0])*omega-648.0*phiC[0]*rdxCp2Sq[1]-2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-720.0*phiC[0]*rdxCp2Sq[0]))/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]-138.5640646055102*rdxCp2[0]*rho[0])*omega*volFac+(((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]+75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1]-216.0*phiC[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiUx[1]+54.0*rdxCp2[0]*phiLy[1]-738.0*rdxCp2[0]*phiC[1]+((-77.94228634059945*phiUy[0])+155.8845726811989*phiUx[0]-77.94228634059945*phiLy[0]+277.1281292110203*bcVals[0])*rdxCp2[0])*rdxCp2[1]-240.0*rdxCp2Sq[0]*phiC[1]+277.1281292110203*bcVals[0]*rdxCp2Sq[0])*omega+216.0*phiC[1]*rdxCp2Sq[1]+738.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+240.0*rdxCp2Sq[0]*phiC[1])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((277.1281292110203*rdxCp2[0]*rho[3]+((-1472.0*rdxCp2[1])-1680.0*rdxCp2[0])*rho[2])*omega*volFac+((-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3])+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-828.0*rdxCp2[0]*rdxCp2[1])-720.0*rdxCp2Sq[0])*phiUx[2]+(644.0*rdxCp2Sq[1]+735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[2]+(796.7433714816835*phiLy[0]-796.7433714816835*phiUy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiLy[0]-909.3266739736605*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[2]))/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[0]*rho[2])*omega*volFac+(((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[3]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]+1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiUx[2]+181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiLy[0]-225.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[3])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[1]-144.0*rho[0]*rdxCp2[1]-1600.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(77.94228634059945*rdxCp2Sq[1]+866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-77.94228634059945*rdxCp2Sq[1])-866.0254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-81.0*phiUy[0])-81.0*phiLy[0]+162.0*phiC[0])*rdxCp2Sq[1]+(420.8883462392369*rdxCp2[0]*phiUy[1]+420.8883462392369*rdxCp2[0]*phiLy[1]+93.53074360871933*rdxCp2[0]*phiLx[1]-864.0*rdxCp2[0]*bcVals[1]+((-900.0*phiUy[0])-900.0*phiLy[0]+54.0*phiLx[0]+2178.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*phiLx[1]-3120.0*rdxCp2Sq[0]*bcVals[1]+(2580.0*phiC[0]-1020.0*phiLx[0])*rdxCp2Sq[0])*omega-162.0*phiC[0]*rdxCp2Sq[1]-2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-2580.0*phiC[0]*rdxCp2Sq[0]))/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[1] = (((144.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[1]-277.1281292110203*rdxCp2[0]*rho[0])*omega*volFac+(((-77.94228634059945*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(77.94228634059945*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*phiLy[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(81.0*phiUy[1]+81.0*phiLy[1]-162.0*phiC[1])*rdxCp2Sq[1]+(189.0*rdxCp2[0]*phiUy[1]+189.0*rdxCp2[0]*phiLy[1]-360.0*rdxCp2[0]*phiLx[1]-2178.0*rdxCp2[0]*phiC[1]+1247.076581449591*rdxCp2[0]*bcVals[1]+((-155.8845726811989*phiUy[0])-155.8845726811989*phiLy[0]-311.7691453623978*phiLx[0])*rdxCp2[0])*rdxCp2[1]-660.0*rdxCp2Sq[0]*phiLx[1]-2580.0*rdxCp2Sq[0]*phiC[1]+1247.076581449591*rdxCp2Sq[0]*bcVals[1]-623.5382907247956*phiLx[0]*rdxCp2Sq[0])*omega+162.0*phiC[1]*rdxCp2Sq[1]+2178.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+2580.0*rdxCp2Sq[0]*phiC[1])/(162.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((748.2459488697547*rdxCp2[0]*rho[3]+((-368.0*rdxCp2[1])-1600.0*rdxCp2[0])*rho[2])*omega*volFac+((-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUy[3])-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0])*phiLx[3]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(161.0*rdxCp2Sq[1]+700.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(138.0*rdxCp2[0]*rdxCp2[1]-1020.0*rdxCp2Sq[0])*phiLx[2]+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[2]+(199.1858428704209*phiLy[0]-199.1858428704209*phiUy[0])*rdxCp2Sq[1]+(405.0*rdxCp2[0]*phiUy[1]-405.0*rdxCp2[0]*phiLy[1]+(866.0254037844386*phiLy[0]-866.0254037844386*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[2]))/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 
  phiC[3] = (((368.0*rdxCp2[1]+336.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[0]*rho[2])*omega*volFac+(((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-161.0*rdxCp2Sq[1])-147.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-920.0*rdxCp2[0]*rdxCp2[1])-660.0*rdxCp2Sq[0])*phiLx[3]+((-1058.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-2580.0*rdxCp2Sq[0])*phiC[3]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[2]+121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-796.7433714816835*rdxCp2[0]*rdxCp2[1])-623.5382907247956*rdxCp2Sq[0])*phiLx[2]+(199.1858428704209*phiUy[1]-199.1858428704209*phiLy[1])*rdxCp2Sq[1]+(181.8653347947321*rdxCp2[0]*phiUy[1]-181.8653347947321*rdxCp2[0]*phiLy[1]+(150.0*phiLy[0]-150.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0])*phiC[3])/(1058.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+2580.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[0]*rho[1]+576.0*rho[0]*rdxCp2[1]+1680.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-311.7691453623978*rdxCp2Sq[1])-909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2Sq[1]+909.3266739736605*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(324.0*phiUy[0]+324.0*phiLy[0]-648.0*phiC[0])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]+155.8845726811989*rdxCp2[0]*phiLy[1]+311.7691453623978*rdxCp2[0]*phiLx[1]+576.0*rdxCp2[0]*bcVals[1]+(945.0*phiUy[0]+945.0*phiLy[0]+324.0*phiLx[0]-2214.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0]*phiLx[1]+2080.0*rdxCp2Sq[0]*bcVals[1]+(720.0*phiLx[0]-720.0*phiC[0])*rdxCp2Sq[0])*omega+648.0*phiC[0]*rdxCp2Sq[1]+2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+720.0*phiC[0]*rdxCp2Sq[0])/(648.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[1] = (((192.0*rdxCp2[1]+96.0*rdxCp2[0])*rho[1]+138.5640646055102*rdxCp2[0]*rho[0])*omega*volFac+(((-103.9230484541326*rdxCp2Sq[1])-51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2Sq[1]+51.96152422706631*rdxCp2[0]*rdxCp2[1])*phiLy[3]-75.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+75.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(108.0*phiUy[1]+108.0*phiLy[1]-216.0*phiC[1])*rdxCp2Sq[1]+(54.0*rdxCp2[0]*phiUy[1]+54.0*rdxCp2[0]*phiLy[1]-150.0*rdxCp2[0]*phiLx[1]-738.0*rdxCp2[0]*phiC[1]+277.1281292110203*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiUy[0]+77.94228634059945*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1]-240.0*rdxCp2Sq[0]*phiC[1]+277.1281292110203*rdxCp2Sq[0]*bcVals[1])*omega+216.0*phiC[1]*rdxCp2Sq[1]+738.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+240.0*rdxCp2Sq[0]*phiC[1])/(216.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+240.0*rdxCp2Sq[0]); 
  phiC[2] = ((277.1281292110203*rdxCp2[0]*rho[3]+(1472.0*rdxCp2[1]+1680.0*rdxCp2[0])*rho[2])*omega*volFac+((-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUy[3])-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(796.7433714816835*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiLx[3]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-644.0*rdxCp2Sq[1])-735.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(828.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiLx[2]+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[2]+(796.7433714816835*phiUy[0]-796.7433714816835*phiLy[0])*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUy[1]-150.0*rdxCp2[0]*phiLy[1]+(909.3266739736605*phiUy[0]-909.3266739736605*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[2])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 
  phiC[3] = (((1472.0*rdxCp2[1]+288.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[0]*rho[2])*omega*volFac+(((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-644.0*rdxCp2Sq[1])-126.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-4232.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-720.0*rdxCp2Sq[0])*phiC[3]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiUy[2]-181.8653347947321*rdxCp2[0]*rdxCp2[1]*phiLy[2]-1195.115057222525*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(796.7433714816835*phiUy[1]-796.7433714816835*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(225.0*phiUy[0]-225.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0])*phiC[3])/(4232.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+720.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((748.2459488697547*rdxCp2[1]*rho[2]+1600.0*rho[0]*rdxCp2[1]+144.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(3120.0*rdxCp2Sq[1]+864.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(1020.0*phiUy[0]-2580.0*phiC[0])*rdxCp2Sq[1]+((-866.0254037844386*rdxCp2[0]*phiUx[1])+866.0254037844386*rdxCp2[0]*phiLx[1]+((-54.0*phiUy[0])+900.0*phiUx[0]+900.0*phiLx[0]-2178.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-77.94228634059945*rdxCp2Sq[0]*phiUx[1]+77.94228634059945*rdxCp2Sq[0]*phiLx[1]+(81.0*phiUx[0]+81.0*phiLx[0]-162.0*phiC[0])*rdxCp2Sq[0])*omega+2580.0*phiC[0]*rdxCp2Sq[1]+2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+162.0*phiC[0]*rdxCp2Sq[0])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = ((748.2459488697547*rdxCp2[1]*rho[3]+(1600.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[1])*omega*volFac+((239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiUy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(1020.0*phiUy[1]-2580.0*phiC[1])*rdxCp2Sq[1]+((-138.0*rdxCp2[0]*phiUy[1])-700.0*rdxCp2[0]*phiUx[1]-700.0*rdxCp2[0]*phiLx[1]-5566.0*rdxCp2[0]*phiC[1]+(866.0254037844386*phiUx[0]-866.0254037844386*phiLx[0])*rdxCp2[0])*rdxCp2[1]-161.0*rdxCp2Sq[0]*phiUx[1]-161.0*rdxCp2Sq[0]*phiLx[1]-1058.0*rdxCp2Sq[0]*phiC[1]+(199.1858428704209*phiUx[0]-199.1858428704209*phiLx[0])*rdxCp2Sq[0])*omega+2580.0*phiC[1]*rdxCp2Sq[1]+5566.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+1058.0*rdxCp2Sq[0]*phiC[1])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]+277.1281292110203*rho[0]*rdxCp2[1])*omega*volFac+(((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]+((-2580.0*rdxCp2Sq[1])-2178.0*rdxCp2[0]*rdxCp2[1]-162.0*rdxCp2Sq[0])*phiC[2]+((-1247.076581449591*rdxCp2Sq[1])-1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[2]+623.5382907247956*phiUy[0]*rdxCp2Sq[1]+((-150.0*rdxCp2[0]*phiUx[1])+150.0*rdxCp2[0]*phiLx[1]+(311.7691453623978*phiUy[0]+155.8845726811989*phiUx[0]+155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0])*phiC[2])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]+277.1281292110203*rdxCp2[1]*rho[1])*omega*volFac+(((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+((-2580.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-1058.0*rdxCp2Sq[0])*phiC[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]+623.5382907247956*phiUy[1]*rdxCp2Sq[1]+(796.7433714816835*rdxCp2[0]*phiUy[1]-121.2435565298214*rdxCp2[0]*phiUx[1]-121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiUx[0]-150.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0])*phiC[3])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[2]-1680.0*rho[0]*rdxCp2[1]-576.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(720.0*phiC[0]-720.0*phiUy[0])*rdxCp2Sq[1]+(909.3266739736605*rdxCp2[0]*phiUx[1]-909.3266739736605*rdxCp2[0]*phiLx[1]+((-324.0*phiUy[0])-945.0*phiUx[0]-945.0*phiLx[0]+2214.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+311.7691453623978*rdxCp2Sq[0]*phiUx[1]-311.7691453623978*rdxCp2Sq[0]*phiLx[1]+((-324.0*phiUx[0])-324.0*phiLx[0]+648.0*phiC[0])*rdxCp2Sq[0])*omega-720.0*phiC[0]*rdxCp2Sq[1]-2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-648.0*phiC[0]*rdxCp2Sq[0]))/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((277.1281292110203*rdxCp2[1]*rho[3]+((-1680.0*rdxCp2[1])-1472.0*rdxCp2[0])*rho[1])*omega*volFac+((692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiUy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(720.0*phiC[1]-720.0*phiUy[1])*rdxCp2Sq[1]+((-828.0*rdxCp2[0]*phiUy[1])+735.0*rdxCp2[0]*phiUx[1]+735.0*rdxCp2[0]*phiLx[1]+5658.0*rdxCp2[0]*phiC[1]+(909.3266739736605*phiLx[0]-909.3266739736605*phiUx[0])*rdxCp2[0])*rdxCp2[1]+644.0*rdxCp2Sq[0]*phiUx[1]+644.0*rdxCp2Sq[0]*phiLx[1]+4232.0*rdxCp2Sq[0]*phiC[1]+(796.7433714816835*phiLx[0]-796.7433714816835*phiUx[0])*rdxCp2Sq[0])*omega-720.0*phiC[1]*rdxCp2Sq[1]-5658.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-4232.0*rdxCp2Sq[0]*phiC[1]))/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]-138.5640646055102*rho[0]*rdxCp2[1])*omega*volFac+(((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]-150.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+((-240.0*rdxCp2Sq[1])-738.0*rdxCp2[0]*rdxCp2[1]-216.0*rdxCp2Sq[0])*phiC[2]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(75.0*rdxCp2[0]*phiUx[1]-75.0*rdxCp2[0]*phiLx[1]+(155.8845726811989*phiUy[0]-77.94228634059945*phiUx[0]-77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0])*phiC[2])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]-415.6921938165305*rdxCp2[1]*rho[1])*omega*volFac+((-1150.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+((-720.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-4232.0*rdxCp2Sq[0])*phiC[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+(1195.115057222525*rdxCp2[0]*phiUy[1]+181.8653347947321*rdxCp2[0]*phiUx[1]+181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiLx[0]-225.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0])*phiC[3])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[2]-1600.0*rho[0]*rdxCp2[1]-144.0*rdxCp2[0]*rho[0])*omega*volFac+((-405.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+405.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-3120.0*rdxCp2Sq[1])-864.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(93.53074360871933*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[2]+420.8883462392369*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2580.0*phiC[0]-1020.0*phiLy[0])*rdxCp2Sq[1]+(866.0254037844386*rdxCp2[0]*phiUx[1]-866.0254037844386*rdxCp2[0]*phiLx[1]+((-900.0*phiUx[0])+54.0*phiLy[0]-900.0*phiLx[0]+2178.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0]*phiUx[1]-77.94228634059945*rdxCp2Sq[0]*phiLx[1]+((-81.0*phiUx[0])-81.0*phiLx[0]+162.0*phiC[0])*rdxCp2Sq[0])*omega-2580.0*phiC[0]*rdxCp2Sq[1]-2178.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-162.0*phiC[0]*rdxCp2Sq[0]))/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((748.2459488697547*rdxCp2[1]*rho[3]+((-1600.0*rdxCp2[1])-368.0*rdxCp2[0])*rho[1])*omega*volFac+((-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiUx[3])+(239.023011444505*rdxCp2[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[1])*phiLy[3]-327.3576026305176*rdxCp2[0]*rdxCp2[1]*phiLx[3]+405.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-405.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(2580.0*phiC[1]-1020.0*phiLy[1])*rdxCp2Sq[1]+(700.0*rdxCp2[0]*phiUx[1]+138.0*rdxCp2[0]*phiLy[1]+700.0*rdxCp2[0]*phiLx[1]+5566.0*rdxCp2[0]*phiC[1]+(866.0254037844386*phiLx[0]-866.0254037844386*phiUx[0])*rdxCp2[0])*rdxCp2[1]+161.0*rdxCp2Sq[0]*phiUx[1]+161.0*rdxCp2Sq[0]*phiLx[1]+1058.0*rdxCp2Sq[0]*phiC[1]+(199.1858428704209*phiLx[0]-199.1858428704209*phiUx[0])*rdxCp2Sq[0])*omega-2580.0*phiC[1]*rdxCp2Sq[1]-5566.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-1058.0*rdxCp2Sq[0]*phiC[1]))/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 
  phiC[2] = (((336.0*rdxCp2[1]+144.0*rdxCp2[0])*rho[2]-277.1281292110203*rho[0]*rdxCp2[1])*omega*volFac+(((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-77.94228634059945*rdxCp2Sq[0])*phiUx[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*phiLx[3]+(1247.076581449591*rdxCp2Sq[1]+1247.076581449591*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiUx[2]+((-660.0*rdxCp2Sq[1])-360.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(189.0*rdxCp2[0]*rdxCp2[1]+81.0*rdxCp2Sq[0])*phiLx[2]+((-2580.0*rdxCp2Sq[1])-2178.0*rdxCp2[0]*rdxCp2[1]-162.0*rdxCp2Sq[0])*phiC[2]-623.5382907247956*phiLy[0]*rdxCp2Sq[1]+(150.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiLx[1]+((-155.8845726811989*phiUx[0])-311.7691453623978*phiLy[0]-155.8845726811989*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0])*phiC[2])/(2580.0*rdxCp2Sq[1]+2178.0*rdxCp2[0]*rdxCp2[1]+162.0*rdxCp2Sq[0]); 
  phiC[3] = (((336.0*rdxCp2[1]+368.0*rdxCp2[0])*rho[3]-277.1281292110203*rdxCp2[1]*rho[1])*omega*volFac+(((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiUx[3]+((-660.0*rdxCp2Sq[1])-920.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-147.0*rdxCp2[0]*rdxCp2[1])-161.0*rdxCp2Sq[0])*phiLx[3]+((-2580.0*rdxCp2Sq[1])-5566.0*rdxCp2[0]*rdxCp2[1]-1058.0*rdxCp2Sq[0])*phiC[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+199.1858428704209*rdxCp2Sq[0])*phiUx[2]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-199.1858428704209*rdxCp2Sq[0])*phiLx[2]-623.5382907247956*phiLy[1]*rdxCp2Sq[1]+(121.2435565298214*rdxCp2[0]*phiUx[1]-796.7433714816835*rdxCp2[0]*phiLy[1]+121.2435565298214*rdxCp2[0]*phiLx[1]+(150.0*phiLx[0]-150.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0])*phiC[3])/(2580.0*rdxCp2Sq[1]+5566.0*rdxCp2[0]*rdxCp2[1]+1058.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((277.1281292110203*rdxCp2[1]*rho[2]+1680.0*rho[0]*rdxCp2[1]+576.0*rdxCp2[0]*rho[0])*omega*volFac+((-150.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+150.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(2080.0*rdxCp2Sq[1]+576.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(692.8203230275509*rdxCp2Sq[1]+311.7691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[2]+155.8845726811989*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(720.0*phiLy[0]-720.0*phiC[0])*rdxCp2Sq[1]+((-909.3266739736605*rdxCp2[0]*phiUx[1])+909.3266739736605*rdxCp2[0]*phiLx[1]+(945.0*phiUx[0]+324.0*phiLy[0]+945.0*phiLx[0]-2214.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-311.7691453623978*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLx[1]+(324.0*phiUx[0]+324.0*phiLx[0]-648.0*phiC[0])*rdxCp2Sq[0])*omega+720.0*phiC[0]*rdxCp2Sq[1]+2214.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+648.0*phiC[0]*rdxCp2Sq[0])/(720.0*rdxCp2Sq[1]+2214.0*rdxCp2[0]*rdxCp2[1]+648.0*rdxCp2Sq[0]); 
  phiC[1] = ((277.1281292110203*rdxCp2[1]*rho[3]+(1680.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[1])*omega*volFac+((-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiUx[3])+(692.8203230275509*rdxCp2Sq[1]+796.7433714816835*rdxCp2[0]*rdxCp2[1])*phiLy[3]-121.2435565298214*rdxCp2[0]*rdxCp2[1]*phiLx[3]+150.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(720.0*phiLy[1]-720.0*phiC[1])*rdxCp2Sq[1]+((-735.0*rdxCp2[0]*phiUx[1])+828.0*rdxCp2[0]*phiLy[1]-735.0*rdxCp2[0]*phiLx[1]-5658.0*rdxCp2[0]*phiC[1]+(909.3266739736605*phiUx[0]-909.3266739736605*phiLx[0])*rdxCp2[0])*rdxCp2[1]-644.0*rdxCp2Sq[0]*phiUx[1]-644.0*rdxCp2Sq[0]*phiLx[1]-4232.0*rdxCp2Sq[0]*phiC[1]+(796.7433714816835*phiUx[0]-796.7433714816835*phiLx[0])*rdxCp2Sq[0])*omega+720.0*phiC[1]*rdxCp2Sq[1]+5658.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+4232.0*rdxCp2Sq[0]*phiC[1])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 
  phiC[2] = (((96.0*rdxCp2[1]+192.0*rdxCp2[0])*rho[2]+138.5640646055102*rho[0]*rdxCp2[1])*omega*volFac+(((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiLx[3]+(277.1281292110203*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiUx[2]-150.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(54.0*rdxCp2[0]*rdxCp2[1]+108.0*rdxCp2Sq[0])*phiLx[2]+((-240.0*rdxCp2Sq[1])-738.0*rdxCp2[0]*rdxCp2[1]-216.0*rdxCp2Sq[0])*phiC[2]+((-75.0*rdxCp2[0]*phiUx[1])+75.0*rdxCp2[0]*phiLx[1]+(77.94228634059945*phiUx[0]-155.8845726811989*phiLy[0]+77.94228634059945*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0])*phiC[2])/(240.0*rdxCp2Sq[1]+738.0*rdxCp2[0]*rdxCp2[1]+216.0*rdxCp2Sq[0]); 
  phiC[3] = (((288.0*rdxCp2[1]+1472.0*rdxCp2[0])*rho[3]+415.6921938165305*rdxCp2[1]*rho[1])*omega*volFac+(((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiUx[3]-1150.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-126.0*rdxCp2[0]*rdxCp2[1])-644.0*rdxCp2Sq[0])*phiLx[3]+((-720.0*rdxCp2Sq[1])-5658.0*rdxCp2[0]*rdxCp2[1]-4232.0*rdxCp2Sq[0])*phiC[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+796.7433714816835*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-796.7433714816835*rdxCp2Sq[0])*phiLx[2]+((-181.8653347947321*rdxCp2[0]*phiUx[1])-1195.115057222525*rdxCp2[0]*phiLy[1]-181.8653347947321*rdxCp2[0]*phiLx[1]+(225.0*phiUx[0]-225.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0])*phiC[3])/(720.0*rdxCp2Sq[1]+5658.0*rdxCp2[0]*rdxCp2[1]+4232.0*rdxCp2Sq[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(328950.0*phiUy[0]-832050.0*phiC[0])*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiUx[1]+(2715915.0*phiUy[0]-193500.0*phiUx[0]-8611395.0*phiC[0]+3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiUx[1]+((-193500.0*phiUy[0])+2715915.0*phiUx[0]-8611395.0*phiC[0]+9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiUx[1]+(328950.0*phiUx[0]-832050.0*phiC[0]+1006200.0*bcVals[0])*rdxCp2R3[0])*omega+832050.0*phiC[0]*rdxCp2R3[1]+8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+832050.0*phiC[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(328950.0*phiUy[1]-832050.0*phiC[1])*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiUx[1]-8611395.0*rdxCp2[0]*phiC[1]+(567939.4598018348*phiUy[0]+1117172.770881926*phiUx[0]-4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiUx[1]-8611395.0*rdxCp2Sq[0]*phiC[1]+((-33515.18312645776*phiUy[0])+1919718.512568964*phiUx[0]-4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiUx[1]-832050.0*rdxCp2R3[0]*phiC[1]+(201091.0987587466*phiUx[0]-402182.1975174932*bcVals[0])*rdxCp2R3[0])*omega+832050.0*phiC[1]*rdxCp2R3[1]+8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+832050.0*rdxCp2R3[0]*phiC[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiUx[1]+(1919718.512568964*phiUy[0]-33515.18312645776*phiUx[0]+536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiUx[1]+(1117172.770881926*phiUy[0]+567939.4598018348*phiUx[0]+1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiUx[1]+(44400.0*phiUy[0]+64500.0*phiUx[0]-258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiUx[1]+(64500.0*phiUy[0]+44400.0*phiUx[0]+40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiC[0]-12960.0*phiUy[0])*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiUx[1]+((-456408.0*phiUy[0])+11340.0*phiUx[0]+535788.0*phiC[0]-181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiUx[1]+((-1842246.0*phiUy[0])-532953.0*phiUx[0]+3688425.0*phiC[0]-2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiUx[1]+((-928800.0*phiUy[0])-2563956.0*phiUx[0]+7679628.0*phiC[0]-8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiUx[1]+((-1052640.0*phiUx[0])+2662560.0*phiC[0]-3219840.0*bcVals[0])*rdxCp2R4[0])*omega-12960.0*phiC[0]*rdxCp2R4[1]-535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-2662560.0*phiC[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+((12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiC[1]-12960.0*phiUy[1])*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiUx[1]+535788.0*rdxCp2[0]*phiC[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiUx[1]+3688425.0*rdxCp2Sq[0]*phiC[1]+((-301792.5327108011*phiUy[0])-659547.6270141526*phiUx[0]+1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiUx[1]+7679628.0*rdxCp2R3[0]*phiC[1]+((-160872.8790069972*phiUy[0])-1745283.675738703*phiUx[0]+3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiUx[1]+2662560.0*rdxCp2R4[0]*phiC[1]+(1286983.032055978*bcVals[0]-643491.516027989*phiUx[0])*rdxCp2R4[0])*omega-12960.0*phiC[1]*rdxCp2R4[1]-535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-2662560.0*rdxCp2R4[0]*phiC[1]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiUx[1]+(98207.28078915528*phiUy[0]+14029.6115413079*phiUx[0]-224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiUx[1]+(2061183.762277152*phiUy[0]-814886.6036909672*phiUx[0]-2492594.31717237*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiUx[1]+(6703036.625291553*phiUy[0]-3407636.758811008*phiUx[0]-6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiUx[1]+(332100.0*bcVals[0]-166050.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiUx[1]+(387000.0*phiUy[0]-266400.0*phiUx[0]-241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(2662560.0*phiC[0]-1052640.0*phiUy[0])*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiUx[1]+((-2563956.0*phiUy[0])-928800.0*phiUx[0]+7679628.0*phiC[0]+1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiUx[1]+((-532953.0*phiUy[0])-1842246.0*phiUx[0]+3688425.0*phiC[0]+5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiUx[1]+(11340.0*phiUy[0]-456408.0*phiUx[0]+535788.0*phiC[0]+1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiUx[1]+((-12960.0*phiUx[0])+12960.0*phiC[0]+37440.0*bcVals[0])*rdxCp2R4[0])*omega-2662560.0*phiC[0]*rdxCp2R4[1]-7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-12960.0*phiC[0]*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(5263200.0*phiUy[1]-1.33128e+7*phiC[1])*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiUx[1]-3.839814e+7*rdxCp2[0]*phiC[1]+((-3407636.758811008*phiUy[0])+6703036.625291553*phiUx[0]+1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiUx[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]+((-814886.6036909672*phiUy[0])+2061183.762277152*phiUx[0]+1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiUx[1]-2678940.0*rdxCp2R3[0]*phiC[1]+(14029.6115413079*phiUy[0]+98207.28078915528*phiUx[0]+2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]+74824.59488697546*bcVals[0]*rdxCp2R4[0])*omega+1.33128e+7*phiC[1]*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiUx[1]+((-1745283.675738703*phiUy[0])-160872.8790069972*phiUx[0]+285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiUx[1]+((-659547.6270141526*phiUy[0])-301792.5327108011*phiUx[0]+871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiUx[1]+((-65471.52052610354*phiUy[0])-65471.52052610354*phiUx[0]+201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiUx[1]+((-266400.0*phiUy[0])+387000.0*phiUx[0]+688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+(580800.0*bcVals[0]-166050.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiUx[1]+((-27000.0*phiUy[0])-27000.0*phiUx[0]+99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiUy[0]-1440.0*phiC[0])*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiUx[1]+(11799.0*phiUy[0]+1890.0*phiUx[0]-13689.0*phiC[0]-3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiUx[1]+(1890.0*phiUy[0]+11799.0*phiUx[0]-13689.0*phiC[0]-33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiUx[1]+(1440.0*phiUx[0]-1440.0*phiC[0]-4160.0*bcVals[0])*rdxCp2R3[0])*omega+1440.0*phiC[0]*rdxCp2R3[1]+13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+1440.0*phiC[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiC[1]-1440.0*phiUy[1])*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiUx[1]+13689.0*rdxCp2[0]*phiC[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiUx[0]-4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiUx[1]+13689.0*rdxCp2Sq[0]*phiC[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiUx[0]-14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1]-1662.768775266122*bcVals[0]*rdxCp2R3[0])*omega-1440.0*phiC[1]*rdxCp2R3[1]-13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiUx[1]+((-467.6537180435967*phiUy[0])+467.6537180435967*phiUx[0]-831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiUx[1]+((-2727.980021920981*phiUy[0])+2727.980021920981*phiUx[0]-8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiUx[1]+(1125.0*phiUy[0]-1125.0*phiUx[0]-2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiUx[1]+((-1125.0*phiUy[0])+1125.0*phiUx[0]-4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(241309.3185104959*rdxCp2Sq[1]+2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1533436.541464953*rdxCp2Sq[0]*rdxCp2[1]-90490.99444143593*rdxCp2[0]*rdxCp2Sq[1])*phiUx[2]+((-268121.4650116621*rdxCp2R3[1])-2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]+335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(832050.0*phiC[0]-328950.0*phiLy[0])*rdxCp2R3[1]+((-335151.8312645776*rdxCp2[0]*phiUx[1])-1533436.541464953*rdxCp2[0]*phiLy[1]+(193500.0*phiUx[0]-2715915.0*phiLy[0]+8611395.0*phiC[0]-3096000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2176434.423012786*rdxCp2Sq[0]*phiUx[1]+90490.99444143593*rdxCp2Sq[0]*phiLy[1]+((-2715915.0*phiUx[0])+193500.0*phiLy[0]+8611395.0*phiC[0]-9081960.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiUx[1]+((-328950.0*phiUx[0])+832050.0*phiC[0]-1006200.0*bcVals[0])*rdxCp2R3[0])*omega-832050.0*phiC[0]*rdxCp2R3[1]-8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-832050.0*phiC[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(522450.0*rdxCp2[0]*rdxCp2Sq[1]+359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(58050.0*rdxCp2Sq[0]*rdxCp2[1]-493650.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+(832050.0*phiC[1]-328950.0*phiLy[1])*rdxCp2R3[1]+(1290000.0*rdxCp2[0]*phiUx[1]-125505.0*rdxCp2[0]*phiLy[1]+8611395.0*rdxCp2[0]*phiC[1]+((-1117172.770881926*phiUx[0])-567939.4598018348*phiLy[0]+4468691.083527703*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2054550.0*rdxCp2Sq[0]*phiUx[1]+40635.0*rdxCp2Sq[0]*phiLy[1]+8611395.0*rdxCp2Sq[0]*phiC[1]+((-1919718.512568964*phiUx[0])+33515.18312645776*phiLy[0]+4308649.588908338*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiUx[1]+832050.0*rdxCp2R3[0]*phiC[1]+(402182.1975174932*bcVals[0]-201091.0987587466*phiUx[0])*rdxCp2R3[0])*omega-832050.0*phiC[1]*rdxCp2R3[1]-8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-832050.0*rdxCp2R3[0]*phiC[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = (((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+(108360.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiUx[3]+((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-40635.0*rdxCp2[0]*rdxCp2Sq[1])+125505.0*rdxCp2Sq[0]*rdxCp2[1]+328950.0*rdxCp2R3[0])*phiUx[2]+((-212850.0*rdxCp2R3[1])-2054550.0*rdxCp2[0]*rdxCp2Sq[1]-1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]-201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-58050.0*rdxCp2[0]*phiUx[1])-359640.0*rdxCp2[0]*phiLy[1]+(33515.18312645776*phiUx[0]-1919718.512568964*phiLy[0]-536242.9300233242*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(493650.0*rdxCp2Sq[0]*phiUx[1]-522450.0*rdxCp2Sq[0]*phiLy[1]+((-567939.4598018348*phiUx[0])-1117172.770881926*phiLy[0]-1098466.622160182*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+(28890.60747024887*rdxCp2[0]*rdxCp2[1]+29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiUx[3]+((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(78202.0939617348*rdxCp2[0]*rdxCp2Sq[1]+437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]+67030.36625291553*rdxCp2R3[0])*phiUx[2]+((-40789.79651824706*rdxCp2[0]*rdxCp2Sq[1])-74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+(74478.18472546172*rdxCp2[0]*phiUx[1]-437394.7904353685*rdxCp2[0]*phiLy[1]+((-64500.0*phiUx[0])-44400.0*phiLy[0]+258000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(40789.79651824706*rdxCp2Sq[0]*phiUx[1]-78202.0939617348*rdxCp2Sq[0]*phiLy[1]+((-44400.0*phiUx[0])-64500.0*phiLy[0]-40200.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(4988.306325798365*rdxCp2R3[1]+170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]+599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(37440.0*rdxCp2R4[1]+1303392.0*rdxCp2[0]*rdxCp2R3[1]+5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-1870.614872174387*rdxCp2[0]*rdxCp2R3[1])+108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]+454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(12470.76581449591*rdxCp2R4[1]+439178.8027671645*rdxCp2[0]*rdxCp2R3[1]+1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]+893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(12960.0*phiLy[0]-12960.0*phiC[0])*rdxCp2R4[1]+(19641.45615783106*rdxCp2[0]*phiUx[1]+176773.1054204795*rdxCp2[0]*phiLy[1]+((-11340.0*phiUx[0])+456408.0*phiLy[0]-535788.0*phiC[0]+181440.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-386469.0325912283*rdxCp2Sq[0]*phiUx[1])+814839.8383191626*rdxCp2Sq[0]*phiLy[1]+(532953.0*phiUx[0]+1842246.0*phiLy[0]-3688425.0*phiC[0]+2626452.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2064285.865273508*rdxCp2R3[0]*phiUx[1])+434356.7733188925*rdxCp2R3[0]*phiLy[1]+(2563956.0*phiUx[0]+928800.0*phiLy[0]-7679628.0*phiC[0]+8373744.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiUx[1]+(1052640.0*phiUx[0]-2662560.0*phiC[0]+3219840.0*bcVals[0])*rdxCp2R4[0])*omega+12960.0*phiC[0]*rdxCp2R4[1]+535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+2662560.0*phiC[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(29520.0*rdxCp2[0]*rdxCp2Sq[1]+116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]+285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(10800.0*rdxCp2[0]*rdxCp2R3[1]+66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]+106560.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(63000.0*rdxCp2[0]*rdxCp2R3[1]+290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(12960.0*phiLy[1]-12960.0*phiC[1])*rdxCp2R4[1]+((-75600.0*rdxCp2[0]*phiUx[1])+157788.0*rdxCp2[0]*phiLy[1]-535788.0*rdxCp2[0]*phiC[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-261886.0821044141*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-727155.0*rdxCp2Sq[0]*phiUx[1])+465750.0*rdxCp2Sq[0]*phiLy[1]-3688425.0*rdxCp2Sq[0]*phiC[1]+(659547.6270141526*phiUx[0]+301792.5327108011*phiLy[0]-1922680.319449907*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1862820.0*rdxCp2R3[0]*phiUx[1])+195048.0*rdxCp2R3[0]*phiLy[1]-7679628.0*rdxCp2R3[0]*phiC[1]+(1745283.675738703*phiUx[0]+160872.8790069972*phiLy[0]-3812313.1094914*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiUx[1]-2662560.0*rdxCp2R4[0]*phiC[1]+(643491.516027989*phiUx[0]-1286983.032055978*bcVals[0])*rdxCp2R4[0])*omega+12960.0*phiC[1]*rdxCp2R4[1]+535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+2662560.0*rdxCp2R4[0]*phiC[1])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = (((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+(25920.0*rdxCp2R3[1]+1006560.0*rdxCp2[0]*rdxCp2Sq[1]+5652000.0*rdxCp2Sq[0]*rdxCp2[1]+8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiUx[3]+(210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(74824.59488697546*rdxCp2R4[1]+2731097.713374604*rdxCp2[0]*rdxCp2R3[1]+1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-9720.0*rdxCp2[0]*rdxCp2R3[1])-63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1408860.0*rdxCp2R3[0]*rdxCp2[1]+5263200.0*rdxCp2R4[0])*phiUx[2]+((-94500.0*rdxCp2[0]*rdxCp2R3[1])-1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]+(24300.0*rdxCp2[0]*phiUx[1]+218700.0*rdxCp2[0]*phiLy[1]+((-14029.6115413079*phiUx[0])-98207.28078915528*phiLy[0]+224473.7846609264*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((814886.6036909672*phiUx[0]-2061183.762277152*phiLy[0]+2492594.31717237*bcVals[0])*rdxCp2Sq[0]-664200.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-2961900.0*rdxCp2R3[0]*phiUx[1])-3134700.0*rdxCp2R3[0]*phiLy[1]+(3407636.758811008*phiUx[0]-6703036.625291553*phiLy[0]+6590799.732961089*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2])/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+(17874.76433411081*rdxCp2[0]*rdxCp2Sq[1]+201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]+476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiUx[3]+((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+(99600.0*rdxCp2[0]*rdxCp2R3[1]+580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(18706.14872174387*rdxCp2[0]*rdxCp2R3[1]+543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]+2016730.678300898*rdxCp2R3[0]*rdxCp2[1]+1072485.860046648*rdxCp2R4[0])*phiUx[2]+(25980.76211353316*rdxCp2[0]*rdxCp2R3[1]-372390.9236273086*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-31176.91453623978*rdxCp2[0]*phiUx[1])-155884.5726811989*rdxCp2[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-108000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-175759.8556980518*rdxCp2Sq[0]*phiUx[1])-687061.2540923841*rdxCp2Sq[0]*phiLy[1]+(166050.0*phiUx[0]-332100.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-244738.7791094823*rdxCp2R3[0]*phiUx[1])-469212.5637704087*rdxCp2R3[0]*phiLy[1]+(266400.0*phiUx[0]-387000.0*phiLy[0]+241200.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-772189.8192335867*rdxCp2R3[1])-1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]-429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+((417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-434356.7733188925*rdxCp2[0]*rdxCp2R3[1])-814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]-176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(857988.6880373188*rdxCp2R4[1]+2064285.865273508*rdxCp2[0]*rdxCp2R3[1]+386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]-19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1052640.0*phiLy[0]-2662560.0*phiC[0])*rdxCp2R4[1]+((-893738.2167055405*rdxCp2[0]*phiUx[1])-454351.5678414679*rdxCp2[0]*phiLy[1]+(928800.0*phiUx[0]+2563956.0*phiLy[0]-7679628.0*phiC[0]-1651200.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-1772702.040022519*rdxCp2Sq[0]*phiUx[1])-108651.5471587956*rdxCp2Sq[0]*phiLy[1]+(1842246.0*phiUx[0]+532953.0*phiLy[0]-3688425.0*phiC[0]-5004704.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-439178.8027671645*rdxCp2R3[0]*phiUx[1])+1870.614872174387*rdxCp2R3[0]*phiLy[1]+(456408.0*phiUx[0]-11340.0*phiLy[0]-535788.0*phiC[0]-1303392.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiUx[1]+(12960.0*phiUx[0]-12960.0*phiC[0]-37440.0*bcVals[0])*rdxCp2R4[0])*omega+2662560.0*phiC[0]*rdxCp2R4[1]+7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+12960.0*phiC[0]*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2352240.0*rdxCp2[0]*rdxCp2Sq[1])-597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+((210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3134700.0*rdxCp2[0]*rdxCp2R3[1]-218700.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(2961900.0*rdxCp2[0]*rdxCp2R3[1]+664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1.33128e+7*phiC[1]-5263200.0*phiLy[1])*rdxCp2R4[1]+(6450000.0*rdxCp2[0]*phiUx[1]-1408860.0*rdxCp2[0]*phiLy[1]+3.839814e+7*rdxCp2[0]*phiC[1]+((-6703036.625291553*phiUx[0])+3407636.758811008*phiLy[0]-1.191650955607388e+7*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(1983375.0*rdxCp2Sq[0]*phiUx[1]+63990.0*rdxCp2Sq[0]*phiLy[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]+((-2061183.762277152*phiUx[0])+814886.6036909672*phiLy[0]-1.26515919188061e+7*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(94500.0*rdxCp2R3[0]*phiUx[1]+9720.0*rdxCp2R3[0]*phiLy[1]+2678940.0*rdxCp2R3[0]*phiC[1]+((-98207.28078915528*phiUx[0])-14029.6115413079*phiLy[0]-2731097.713374604*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1]-74824.59488697546*bcVals[0]*rdxCp2R4[0])*omega-1.33128e+7*phiC[1]*rdxCp2R4[1]-3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+((-346752.0*rdxCp2R3[1])-928080.0*rdxCp2[0]*rdxCp2Sq[1]-332172.0*rdxCp2Sq[0]*rdxCp2[1]-30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiUx[3]+((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-195048.0*rdxCp2[0]*rdxCp2R3[1])-465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]-157788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiUx[2]+(681120.0*rdxCp2R4[1]+1862820.0*rdxCp2[0]*rdxCp2R3[1]+727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]+75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2]+643491.516027989*phiLy[0]*rdxCp2R4[1]+((-154800.0*rdxCp2[0]*phiUx[1])-106560.0*rdxCp2[0]*phiLy[1]+(160872.8790069972*phiUx[0]+1745283.675738703*phiLy[0]-285996.229345773*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-290400.0*rdxCp2Sq[0]*phiUx[1])-66420.0*rdxCp2Sq[0]*phiLy[1]+(301792.5327108011*phiUx[0]+659547.6270141526*phiLy[0]-871845.0944978701*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-63000.0*rdxCp2R3[0]*phiUx[1])-10800.0*rdxCp2R3[0]*phiLy[1]+(65471.52052610354*phiUx[0]+65471.52052610354*phiLy[0]-201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+((-173343.6448214932*rdxCp2[0]*rdxCp2Sq[1])-89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]-12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(469212.5637704087*rdxCp2[0]*rdxCp2R3[1]+687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]+155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(244738.7791094823*rdxCp2[0]*rdxCp2R3[1]+175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]+31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+(372390.9236273086*rdxCp2[0]*phiUx[1]-2016730.678300898*rdxCp2[0]*phiLy[1]+((-387000.0*phiUx[0])+266400.0*phiLy[0]-688000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((166050.0*phiLy[0]-580800.0*bcVals[0])*rdxCp2Sq[0]-543205.7742697513*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-25980.76211353316*rdxCp2R3[0]*phiUx[1])-18706.14872174387*rdxCp2R3[0]*phiLy[1]+(27000.0*phiUx[0]+27000.0*phiLy[0]-99600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-554.2562584220407*rdxCp2Sq[1])-4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-4160.0*rdxCp2R3[1])-33726.0*rdxCp2[0]*rdxCp2Sq[1]-3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-311.7691453623978*rdxCp2[0]*rdxCp2Sq[1])-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1385.640646055102*rdxCp2R3[1])-11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]-1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(1440.0*phiC[0]-1440.0*phiLy[0])*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiUx[1]+1818.653347947321*rdxCp2[0]*phiLy[1]+((-1890.0*phiUx[0])-11799.0*phiLy[0]+13689.0*phiC[0]+3360.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(11353.59304361399*rdxCp2Sq[0]*phiUx[1]+311.7691453623978*rdxCp2Sq[0]*phiLy[1]+((-11799.0*phiUx[0])-1890.0*phiLy[0]+13689.0*phiC[0]+33726.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiUx[1]+((-1440.0*phiUx[0])+1440.0*phiC[0]+4160.0*bcVals[0])*rdxCp2R3[0])*omega-1440.0*phiC[0]*rdxCp2R3[1]-13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-1440.0*phiC[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-8400.446416709054*rdxCp2[0]*rdxCp2Sq[1])-831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(450.0*rdxCp2[0]*rdxCp2Sq[1]-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-2625.0*rdxCp2[0]*rdxCp2Sq[1])-450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(1440.0*phiLy[1]-1440.0*phiC[1])*rdxCp2R3[1]+((-2625.0*rdxCp2[0]*phiUx[1])+2664.0*rdxCp2[0]*phiLy[1]-13689.0*rdxCp2[0]*phiC[1]+(2727.980021920981*phiUx[0]-2727.980021920981*phiLy[0]+4849.742261192856*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiUx[1])+324.0*rdxCp2Sq[0]*phiLy[1]-13689.0*rdxCp2Sq[0]*phiC[1]+(467.6537180435967*phiUx[0]-467.6537180435967*phiLy[0]+14081.57306553497*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]+1662.768775266122*bcVals[0]*rdxCp2R3[0])*omega+1440.0*phiC[1]*rdxCp2R3[1]+13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+((-576.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiUx[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1662.768775266122*rdxCp2R3[1])-14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]-4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-324.0*rdxCp2[0]*rdxCp2Sq[1])-2664.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiUx[2]+(450.0*rdxCp2[0]*rdxCp2Sq[1]+2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2]+(450.0*rdxCp2[0]*phiUx[1]+450.0*rdxCp2[0]*phiLy[1]+((-467.6537180435967*phiUx[0])+467.6537180435967*phiLy[0]+831.384387633061*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2625.0*rdxCp2Sq[0]*phiUx[1]-450.0*rdxCp2Sq[0]*phiLy[1]+((-2727.980021920981*phiUx[0])+2727.980021920981*phiLy[0]+8400.446416709054*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+((-744.7818472546172*rdxCp2[0]*rdxCp2[1])-1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+((-4150.0*rdxCp2[0]*rdxCp2Sq[1])-2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(779.4228634059946*rdxCp2[0]*rdxCp2Sq[1]+4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(1082.531754730548*rdxCp2Sq[0]*rdxCp2[1]-1082.531754730548*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-1082.531754730548*rdxCp2[0]*phiUx[1])-4546.633369868302*rdxCp2[0]*phiLy[1]+(1125.0*phiUx[0]-1125.0*phiLy[0]+2000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1082.531754730548*rdxCp2Sq[0]*phiUx[1]-779.4228634059946*rdxCp2Sq[0]*phiLy[1]+((-1125.0*phiUx[0])+1125.0*phiLy[0]+4150.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+(2022134.676820512*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[1]-516000.0*rho[0]*rdxCp2Sq[1]-4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiUy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1006200.0*rdxCp2R3[1])-9081960.0*rdxCp2[0]*rdxCp2Sq[1]-3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(832050.0*phiC[0]-328950.0*phiUy[0])*rdxCp2R3[1]+(1533436.541464953*rdxCp2[0]*phiUy[1]+335151.8312645776*rdxCp2[0]*phiLx[1]-3096000.0*rdxCp2[0]*bcVals[1]+((-2715915.0*phiUy[0])+193500.0*phiLx[0]+8611395.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90490.99444143593*rdxCp2Sq[0]*phiUy[1])-2176434.423012786*rdxCp2Sq[0]*phiLx[1]-9081960.0*rdxCp2Sq[0]*bcVals[1]+(193500.0*phiUy[0]-2715915.0*phiLx[0]+8611395.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-268121.4650116621*rdxCp2R3[0]*phiLx[1]-1006200.0*rdxCp2R3[0]*bcVals[1]+(832050.0*phiC[0]-328950.0*phiLx[0])*rdxCp2R3[0])*omega-832050.0*phiC[0]*rdxCp2R3[1]-8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-832050.0*phiC[0]*rdxCp2R3[0]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = (((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(516000.0*rdxCp2Sq[1]+1016400.0*rdxCp2[0]*rdxCp2[1]+108360.0*rdxCp2Sq[0])*rho[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]-89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1098466.622160182*rdxCp2[0]*rdxCp2Sq[1])-536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(328950.0*phiUy[1]-832050.0*phiC[1])*rdxCp2R3[1]+(125505.0*rdxCp2[0]*phiUy[1]-1290000.0*rdxCp2[0]*phiLx[1]-8611395.0*rdxCp2[0]*phiC[1]+4468691.083527703*rdxCp2[0]*bcVals[1]+((-567939.4598018348*phiUy[0])-1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-40635.0*rdxCp2Sq[0]*phiUy[1])-2054550.0*rdxCp2Sq[0]*phiLx[1]-8611395.0*rdxCp2Sq[0]*phiC[1]+4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(33515.18312645776*phiUy[0]-1919718.512568964*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-212850.0*rdxCp2R3[0]*phiLx[1]-832050.0*rdxCp2R3[0]*phiC[1]+402182.1975174932*rdxCp2R3[0]*bcVals[1]-201091.0987587466*phiLx[0]*rdxCp2R3[0])*omega+832050.0*phiC[1]*rdxCp2R3[1]+8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+832050.0*rdxCp2R3[0]*phiC[1])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]+392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]-89373.82167055405*rho[0]*rdxCp2Sq[1]-748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2]+(402182.1975174932*rdxCp2R3[1]+4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]+4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-201091.0987587466*phiUy[0]*rdxCp2R3[1]+(359640.0*rdxCp2[0]*phiUy[1]+58050.0*rdxCp2[0]*phiLx[1]-536242.9300233242*rdxCp2[0]*bcVals[1]+(33515.18312645776*phiLx[0]-1919718.512568964*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(522450.0*rdxCp2Sq[0]*phiUy[1]-493650.0*rdxCp2Sq[0]*phiLx[1]-1098466.622160182*rdxCp2Sq[0]*bcVals[1]+((-1117172.770881926*phiUy[0])-567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+(29791.27389018469*rdxCp2Sq[1]+28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]-48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]+(258000.0*rdxCp2Sq[0]*rdxCp2[1]-40200.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+67030.36625291553*phiUy[1]*rdxCp2R3[1]+(437394.7904353685*rdxCp2[0]*phiUy[1]-74478.18472546172*rdxCp2[0]*phiLx[1]+258000.0*rdxCp2[0]*bcVals[1]+((-44400.0*phiUy[0])-64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(78202.0939617348*rdxCp2Sq[0]*phiUy[1]-40789.79651824706*rdxCp2Sq[0]*phiLx[1]-40200.0*rdxCp2Sq[0]*bcVals[1]+((-64500.0*phiUy[0])-44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-429306.1131640216*rdxCp2[0]*rdxCp2Sq[1])-1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]-772189.8192335867*rdxCp2R3[0])*rho[1]+30240.0*rho[0]*rdxCp2R3[1]+1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiUy[0]-12960.0*phiC[0])*rdxCp2R4[1]+((-176773.1054204795*rdxCp2[0]*phiUy[1])-19641.45615783106*rdxCp2[0]*phiLx[1]+181440.0*rdxCp2[0]*bcVals[1]+(456408.0*phiUy[0]-11340.0*phiLx[0]-535788.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-814839.8383191626*rdxCp2Sq[0]*phiUy[1])+386469.0325912283*rdxCp2Sq[0]*phiLx[1]+2626452.0*rdxCp2Sq[0]*bcVals[1]+(1842246.0*phiUy[0]+532953.0*phiLx[0]-3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-434356.7733188925*rdxCp2R3[0]*phiUy[1])+2064285.865273508*rdxCp2R3[0]*phiLx[1]+8373744.0*rdxCp2R3[0]*bcVals[1]+(928800.0*phiUy[0]+2563956.0*phiLx[0]-7679628.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+857988.6880373188*rdxCp2R4[0]*phiLx[1]+3219840.0*rdxCp2R4[0]*bcVals[1]+(1052640.0*phiLx[0]-2662560.0*phiC[0])*rdxCp2R4[0])*omega+12960.0*phiC[0]*rdxCp2R4[1]+535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+2662560.0*phiC[0]*rdxCp2R4[0])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-30240.0*rdxCp2R3[1])-332172.0*rdxCp2[0]*rdxCp2Sq[1]-928080.0*rdxCp2Sq[0]*rdxCp2[1]-346752.0*rdxCp2R3[0])*rho[1]+159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+((12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(12960.0*phiC[1]-12960.0*phiUy[1])*rdxCp2R4[1]+((-157788.0*rdxCp2[0]*phiUy[1])+75600.0*rdxCp2[0]*phiLx[1]+535788.0*rdxCp2[0]*phiC[1]-261886.0821044141*rdxCp2[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-465750.0*rdxCp2Sq[0]*phiUy[1])+727155.0*rdxCp2Sq[0]*phiLx[1]+3688425.0*rdxCp2Sq[0]*phiC[1]-1922680.319449907*rdxCp2Sq[0]*bcVals[1]+(301792.5327108011*phiUy[0]+659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-195048.0*rdxCp2R3[0]*phiUy[1])+1862820.0*rdxCp2R3[0]*phiLx[1]+7679628.0*rdxCp2R3[0]*phiC[1]-3812313.1094914*rdxCp2R3[0]*bcVals[1]+(160872.8790069972*phiUy[0]+1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+681120.0*rdxCp2R4[0]*phiLx[1]+2662560.0*rdxCp2R4[0]*phiC[1]-1286983.032055978*rdxCp2R4[0]*bcVals[1]+643491.516027989*phiLx[0]*rdxCp2R4[0])*omega-12960.0*phiC[1]*rdxCp2R4[1]-535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-2662560.0*rdxCp2R4[0]*phiC[1]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+((-597780.0*rdxCp2[0]*rdxCp2Sq[1])-2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+37412.29744348773*rho[0]*rdxCp2R3[1]+1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-218700.0*rdxCp2[0]*phiUy[1])-24300.0*rdxCp2[0]*phiLx[1]+224473.7846609264*rdxCp2[0]*bcVals[1]+((-98207.28078915528*phiUy[0])-14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(664200.0*rdxCp2Sq[0]*phiLx[1]+2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiLx[0]-2061183.762277152*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3134700.0*rdxCp2R3[0]*phiUy[1]+2961900.0*rdxCp2R3[0]*phiLx[1]+6590799.732961089*rdxCp2R3[0]*bcVals[1]+(3407636.758811008*phiLx[0]-6703036.625291553*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+((-12470.76581449591*rdxCp2R3[1])-89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]-173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(155884.5726811989*rdxCp2[0]*phiUy[1]+31176.91453623978*rdxCp2[0]*phiLx[1]-108000.0*rdxCp2[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(687061.2540923841*rdxCp2Sq[0]*phiUy[1]+175759.8556980518*rdxCp2Sq[0]*phiLx[1]-332100.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(469212.5637704087*rdxCp2R3[0]*phiUy[1]+244738.7791094823*rdxCp2R3[0]*phiLx[1]+241200.0*rdxCp2R3[0]*bcVals[1]+(266400.0*phiLx[0]-387000.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = (((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(599151.015354226*rdxCp2[0]*rdxCp2Sq[1]+170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[1]+1651200.0*rho[0]*rdxCp2R3[1]+4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(3219840.0*rdxCp2R4[1]+8373744.0*rdxCp2[0]*rdxCp2R3[1]+2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]+181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(1052640.0*phiUy[0]-2662560.0*phiC[0])*rdxCp2R4[1]+(454351.5678414679*rdxCp2[0]*phiUy[1]+893738.2167055405*rdxCp2[0]*phiLx[1]+1651200.0*rdxCp2[0]*bcVals[1]+(2563956.0*phiUy[0]+928800.0*phiLx[0]-7679628.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(108651.5471587956*rdxCp2Sq[0]*phiUy[1]+1772702.040022519*rdxCp2Sq[0]*phiLx[1]+5004704.0*rdxCp2Sq[0]*bcVals[1]+(532953.0*phiUy[0]+1842246.0*phiLx[0]-3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-1870.614872174387*rdxCp2R3[0]*phiUy[1])+439178.8027671645*rdxCp2R3[0]*phiLx[1]+1303392.0*rdxCp2R3[0]*bcVals[1]+((-11340.0*phiUy[0])+456408.0*phiLx[0]-535788.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+12470.76581449591*rdxCp2R4[0]*phiLx[1]+37440.0*rdxCp2R4[0]*bcVals[1]+(12960.0*phiLx[0]-12960.0*phiC[0])*rdxCp2R4[0])*omega+2662560.0*phiC[0]*rdxCp2R4[1]+7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+12960.0*phiC[0]*rdxCp2R4[0])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = (((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(8256000.0*rdxCp2R3[1]+5652000.0*rdxCp2[0]*rdxCp2Sq[1]+1006560.0*rdxCp2Sq[0]*rdxCp2[1]+25920.0*rdxCp2R3[0])*rho[1]+4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(6590799.732961089*rdxCp2[0]*rdxCp2R3[1]+2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]+224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(5263200.0*phiUy[1]-1.33128e+7*phiC[1])*rdxCp2R4[1]+(1408860.0*rdxCp2[0]*phiUy[1]-6450000.0*rdxCp2[0]*phiLx[1]-3.839814e+7*rdxCp2[0]*phiC[1]+1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(3407636.758811008*phiUy[0]-6703036.625291553*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-63990.0*rdxCp2Sq[0]*phiUy[1])-1983375.0*rdxCp2Sq[0]*phiLx[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]+1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(814886.6036909672*phiUy[0]-2061183.762277152*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-9720.0*rdxCp2R3[0]*phiUy[1])-94500.0*rdxCp2R3[0]*phiLx[1]-2678940.0*rdxCp2R3[0]*phiC[1]+2731097.713374604*rdxCp2R3[0]*bcVals[1]+((-14029.6115413079*phiUy[0])-98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]+74824.59488697546*rdxCp2R4[0]*bcVals[1])*omega+1.33128e+7*phiC[1]*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1])/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+(116160.0*rdxCp2[0]*rdxCp2Sq[1]+29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+285996.229345773*rho[0]*rdxCp2R3[1]+704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]+((-1286983.032055978*rdxCp2R4[1])-3812313.1094914*rdxCp2[0]*rdxCp2R3[1]-1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]-261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+643491.516027989*phiUy[0]*rdxCp2R4[1]+(106560.0*rdxCp2[0]*phiUy[1]+154800.0*rdxCp2[0]*phiLx[1]+285996.229345773*rdxCp2[0]*bcVals[1]+(1745283.675738703*phiUy[0]+160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(66420.0*rdxCp2Sq[0]*phiUy[1]+290400.0*rdxCp2Sq[0]*phiLx[1]+871845.0944978701*rdxCp2Sq[0]*bcVals[1]+(659547.6270141526*phiUy[0]+301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(10800.0*rdxCp2R3[0]*phiUy[1]+63000.0*rdxCp2R3[0]*phiLx[1]+201610.7140010173*rdxCp2R3[0]*bcVals[1]+(65471.52052610354*phiUy[0]+65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+(476660.382242955*rdxCp2R3[1]+201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]+17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(241200.0*rdxCp2[0]*rdxCp2R3[1]-332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]-108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1072485.860046648*phiUy[1]*rdxCp2R4[1]+(2016730.678300898*rdxCp2[0]*phiUy[1]-372390.9236273086*rdxCp2[0]*phiLx[1]+688000.0*rdxCp2[0]*bcVals[1]+(266400.0*phiUy[0]-387000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(543205.7742697513*rdxCp2Sq[0]*phiUy[1]+580800.0*rdxCp2Sq[0]*bcVals[1]+166050.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(18706.14872174387*rdxCp2R3[0]*phiUy[1]+25980.76211353316*rdxCp2R3[0]*phiLx[1]+99600.0*rdxCp2R3[0]*bcVals[1]+(27000.0*phiUy[0]+27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4416.729559300637*rdxCp2[0]*rdxCp2[1])-554.2562584220407*rdxCp2Sq[0])*rho[1]-3360.0*rho[0]*rdxCp2Sq[1]-27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]-3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiC[0]-1440.0*phiUy[0])*rdxCp2R3[1]+((-1818.653347947321*rdxCp2[0]*phiUy[1])-1818.653347947321*rdxCp2[0]*phiLx[1]-3360.0*rdxCp2[0]*bcVals[1]+((-11799.0*phiUy[0])-1890.0*phiLx[0]+13689.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2Sq[0]*phiUy[1])-11353.59304361399*rdxCp2Sq[0]*phiLx[1]-33726.0*rdxCp2Sq[0]*bcVals[1]+((-1890.0*phiUy[0])-11799.0*phiLx[0]+13689.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-1385.640646055102*rdxCp2R3[0]*phiLx[1]-4160.0*rdxCp2R3[0]*bcVals[1]+(1440.0*phiC[0]-1440.0*phiLx[0])*rdxCp2R3[0])*omega-1440.0*phiC[0]*rdxCp2R3[1]-13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-1440.0*phiC[0]*rdxCp2R3[0]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-3360.0*rdxCp2Sq[1])-5166.0*rdxCp2[0]*rdxCp2[1]-576.0*rdxCp2Sq[0])*rho[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]-831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1440.0*phiC[1]-1440.0*phiUy[1])*rdxCp2R3[1]+((-2664.0*rdxCp2[0]*phiUy[1])+2625.0*rdxCp2[0]*phiLx[1]+13689.0*rdxCp2[0]*phiC[1]-4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-324.0*rdxCp2Sq[0]*phiUy[1])+450.0*rdxCp2Sq[0]*phiLx[1]+13689.0*rdxCp2Sq[0]*phiC[1]-14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1]-1662.768775266122*rdxCp2R3[0]*bcVals[1])*omega-1440.0*phiC[1]*rdxCp2R3[1]-13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]))/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]-1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]-831.384387633061*rho[0]*rdxCp2Sq[1]-6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-450.0*rdxCp2[0]*phiUy[1])-450.0*rdxCp2[0]*phiLx[1]-831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiUy[0]-467.6537180435967*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(450.0*rdxCp2Sq[0]*phiUy[1]-2625.0*rdxCp2Sq[0]*phiLx[1]-8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiUy[0]-2727.980021920981*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+((-1385.640646055102*rdxCp2Sq[1])-744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]-3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(4546.633369868302*rdxCp2[0]*phiUy[1]+1082.531754730548*rdxCp2[0]*phiLx[1]-2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(779.4228634059946*rdxCp2Sq[0]*phiUy[1]-1082.531754730548*rdxCp2Sq[0]*phiLx[1]-4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiUy[0]-1125.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((1058508.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-241309.3185104959*rdxCp2Sq[1])-2022134.676820512*rdxCp2[0]*rdxCp2[1])*rho[2]+((-2022134.676820512*rdxCp2[0]*rdxCp2[1])-241309.3185104959*rdxCp2Sq[0])*rho[1]+516000.0*rho[0]*rdxCp2Sq[1]+4432360.0*rdxCp2[0]*rho[0]*rdxCp2[1]+516000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((156735.0*rdxCp2Sq[0]*rdxCp2[1]-1332855.0*rdxCp2[0]*rdxCp2Sq[1])*phiLy[3]+(156735.0*rdxCp2[0]*rdxCp2Sq[1]-1332855.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1006200.0*rdxCp2R3[1]+9081960.0*rdxCp2[0]*rdxCp2Sq[1]+3096000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(268121.4650116621*rdxCp2R3[1]+2176434.423012786*rdxCp2[0]*rdxCp2Sq[1]-335151.8312645776*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(90490.99444143593*rdxCp2[0]*rdxCp2Sq[1]-1533436.541464953*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(328950.0*phiLy[0]-832050.0*phiC[0])*rdxCp2R3[1]+((-1533436.541464953*rdxCp2[0]*phiLy[1])-335151.8312645776*rdxCp2[0]*phiLx[1]+3096000.0*rdxCp2[0]*bcVals[1]+(2715915.0*phiLy[0]-193500.0*phiLx[0]-8611395.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(90490.99444143593*rdxCp2Sq[0]*phiLy[1]+2176434.423012786*rdxCp2Sq[0]*phiLx[1]+9081960.0*rdxCp2Sq[0]*bcVals[1]+((-193500.0*phiLy[0])+2715915.0*phiLx[0]-8611395.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+268121.4650116621*rdxCp2R3[0]*phiLx[1]+1006200.0*rdxCp2R3[0]*bcVals[1]+(328950.0*phiLx[0]-832050.0*phiC[0])*rdxCp2R3[0])*omega+832050.0*phiC[0]*rdxCp2R3[1]+8611395.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+832050.0*phiC[0]*rdxCp2R3[0])/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((241309.3185104959*rdxCp2Sq[1]+234013.9205090157*rdxCp2[0]*rdxCp2[1])*rho[3]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-516000.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-108360.0*rdxCp2Sq[0])*rho[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1]+89373.82167055405*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-268121.4650116621*rdxCp2R3[1])+75136.36403233788*rdxCp2[0]*rdxCp2Sq[1]+70381.8845655613*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-603273.2962762396*rdxCp2[0]*rdxCp2Sq[1])-330397.351797801*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1098466.622160182*rdxCp2[0]*rdxCp2Sq[1]+536242.9300233242*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(493650.0*rdxCp2[0]*rdxCp2Sq[1]-58050.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-522450.0*rdxCp2[0]*rdxCp2Sq[1])-359640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(832050.0*phiC[1]-328950.0*phiLy[1])*rdxCp2R3[1]+((-125505.0*rdxCp2[0]*phiLy[1])+1290000.0*rdxCp2[0]*phiLx[1]+8611395.0*rdxCp2[0]*phiC[1]-4468691.083527703*rdxCp2[0]*bcVals[1]+(567939.4598018348*phiLy[0]+1117172.770881926*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(40635.0*rdxCp2Sq[0]*phiLy[1]+2054550.0*rdxCp2Sq[0]*phiLx[1]+8611395.0*rdxCp2Sq[0]*phiC[1]-4308649.588908338*rdxCp2Sq[0]*bcVals[1]+(1919718.512568964*phiLx[0]-33515.18312645776*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1]+212850.0*rdxCp2R3[0]*phiLx[1]+832050.0*rdxCp2R3[0]*phiC[1]-402182.1975174932*rdxCp2R3[0]*bcVals[1]+201091.0987587466*phiLx[0]*rdxCp2R3[0])*omega-832050.0*phiC[1]*rdxCp2R3[1]-8611395.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-832050.0*rdxCp2R3[0]*phiC[1]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((234013.9205090157*rdxCp2[0]*rdxCp2[1]+241309.3185104959*rdxCp2Sq[0])*rho[3]+((-108360.0*rdxCp2Sq[1])-1016400.0*rdxCp2[0]*rdxCp2[1]-516000.0*rdxCp2Sq[0])*rho[2]-392040.0*rdxCp2[0]*rdxCp2[1]*rho[1]+89373.82167055405*rho[0]*rdxCp2Sq[1]+748938.7691927825*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-330397.351797801*rdxCp2[0]*rdxCp2Sq[1])-603273.2962762396*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(70381.8845655613*rdxCp2[0]*rdxCp2Sq[1]+75136.36403233788*rdxCp2Sq[0]*rdxCp2[1]-268121.4650116621*rdxCp2R3[0])*phiLx[3]+((-402182.1975174932*rdxCp2R3[1])-4308649.588908338*rdxCp2[0]*rdxCp2Sq[1]-4468691.083527703*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(212850.0*rdxCp2R3[1]+2054550.0*rdxCp2[0]*rdxCp2Sq[1]+1290000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(40635.0*rdxCp2[0]*rdxCp2Sq[1]-125505.0*rdxCp2Sq[0]*rdxCp2[1]-328950.0*rdxCp2R3[0])*phiLx[2]+(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0])*phiC[2]+201091.0987587466*phiLy[0]*rdxCp2R3[1]+((-359640.0*rdxCp2[0]*phiLy[1])-58050.0*rdxCp2[0]*phiLx[1]+536242.9300233242*rdxCp2[0]*bcVals[1]+(1919718.512568964*phiLy[0]-33515.18312645776*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-522450.0*rdxCp2Sq[0]*phiLy[1])+493650.0*rdxCp2Sq[0]*phiLx[1]+1098466.622160182*rdxCp2Sq[0]*bcVals[1]+(1117172.770881926*phiLy[0]+567939.4598018348*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-832050.0*rdxCp2R3[1])-8611395.0*rdxCp2[0]*rdxCp2Sq[1]-8611395.0*rdxCp2Sq[0]*rdxCp2[1]-832050.0*rdxCp2R3[0])*phiC[2]))/(832050.0*rdxCp2R3[1]+8611395.0*rdxCp2[0]*rdxCp2Sq[1]+8611395.0*rdxCp2Sq[0]*rdxCp2[1]+832050.0*rdxCp2R3[0]); 
  phiC[3] = (((36120.0*rdxCp2Sq[1]+207028.0*rdxCp2[0]*rdxCp2[1]+36120.0*rdxCp2Sq[0])*rho[3]+((-28890.60747024887*rdxCp2[0]*rdxCp2[1])-29791.27389018469*rdxCp2Sq[0])*rho[2]+((-29791.27389018469*rdxCp2Sq[1])-28890.60747024887*rdxCp2[0]*rdxCp2[1])*rho[1]+48400.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-70950.0*rdxCp2R3[1])-498805.0*rdxCp2[0]*rdxCp2Sq[1]-90300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-90300.0*rdxCp2[0]*rdxCp2Sq[1])-498805.0*rdxCp2Sq[0]*rdxCp2[1]-70950.0*rdxCp2R3[0])*phiLx[3]+((-277350.0*rdxCp2R3[1])-2870465.0*rdxCp2[0]*rdxCp2Sq[1]-2870465.0*rdxCp2Sq[0]*rdxCp2[1]-277350.0*rdxCp2R3[0])*phiC[3]+(40200.0*rdxCp2[0]*rdxCp2Sq[1]-258000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(40789.79651824706*rdxCp2[0]*rdxCp2Sq[1]+74478.18472546172*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-78202.0939617348*rdxCp2[0]*rdxCp2Sq[1])-437394.7904353685*rdxCp2Sq[0]*rdxCp2[1]-67030.36625291553*rdxCp2R3[0])*phiLx[2]-67030.36625291553*phiLy[1]*rdxCp2R3[1]+((-437394.7904353685*rdxCp2[0]*phiLy[1])+74478.18472546172*rdxCp2[0]*phiLx[1]-258000.0*rdxCp2[0]*bcVals[1]+(44400.0*phiLy[0]+64500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-78202.0939617348*rdxCp2Sq[0]*phiLy[1])+40789.79651824706*rdxCp2Sq[0]*phiLx[1]+40200.0*rdxCp2Sq[0]*bcVals[1]+(64500.0*phiLy[0]+44400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0])*phiC[3])/(277350.0*rdxCp2R3[1]+2870465.0*rdxCp2[0]*rdxCp2Sq[1]+2870465.0*rdxCp2Sq[0]*rdxCp2[1]+277350.0*rdxCp2R3[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((79704.0*rdxCp2[0]*rdxCp2Sq[1]+313632.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-4988.306325798365*rdxCp2R3[1])-170433.7994647775*rdxCp2[0]*rdxCp2Sq[1]-599151.015354226*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(429306.1131640216*rdxCp2[0]*rdxCp2Sq[1]+1901292.956078046*rdxCp2Sq[0]*rdxCp2[1]+772189.8192335867*rdxCp2R3[0])*rho[1]-30240.0*rho[0]*rdxCp2R3[1]-1057392.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4139904.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1651200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((170100.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+417960.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3240.0*rdxCp2[0]*rdxCp2R3[1]-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-394920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-37440.0*rdxCp2R4[1])-1303392.0*rdxCp2[0]*rdxCp2R3[1]-5004704.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1651200.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-12470.76581449591*rdxCp2R4[1])-439178.8027671645*rdxCp2[0]*rdxCp2R3[1]-1772702.040022519*rdxCp2Sq[0]*rdxCp2Sq[1]-893738.2167055405*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1870.614872174387*rdxCp2[0]*rdxCp2R3[1]-108651.5471587956*rdxCp2Sq[0]*rdxCp2Sq[1]-454351.5678414679*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(12960.0*phiC[0]-12960.0*phiLy[0])*rdxCp2R4[1]+(176773.1054204795*rdxCp2[0]*phiLy[1]+19641.45615783106*rdxCp2[0]*phiLx[1]-181440.0*rdxCp2[0]*bcVals[1]+((-456408.0*phiLy[0])+11340.0*phiLx[0]+535788.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(814839.8383191626*rdxCp2Sq[0]*phiLy[1]-386469.0325912283*rdxCp2Sq[0]*phiLx[1]-2626452.0*rdxCp2Sq[0]*bcVals[1]+((-1842246.0*phiLy[0])-532953.0*phiLx[0]+3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(434356.7733188925*rdxCp2R3[0]*phiLy[1]-2064285.865273508*rdxCp2R3[0]*phiLx[1]-8373744.0*rdxCp2R3[0]*bcVals[1]+((-928800.0*phiLy[0])-2563956.0*phiLx[0]+7679628.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-857988.6880373188*rdxCp2R4[0]*phiLx[1]-3219840.0*rdxCp2R4[0]*bcVals[1]+(2662560.0*phiC[0]-1052640.0*phiLx[0])*rdxCp2R4[0])*omega-12960.0*phiC[0]*rdxCp2R4[1]-535788.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-7679628.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-2662560.0*phiC[0]*rdxCp2R4[0]))/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[1] = (((4988.306325798365*rdxCp2R3[1]+35791.09788760327*rdxCp2[0]*rdxCp2Sq[1]+69337.45792859727*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-29520.0*rdxCp2[0]*rdxCp2Sq[1])-116160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(30240.0*rdxCp2R3[1]+332172.0*rdxCp2[0]*rdxCp2Sq[1]+928080.0*rdxCp2Sq[0]*rdxCp2[1]+346752.0*rdxCp2R3[0])*rho[1]-159002.2641348229*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-704182.5763252026*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-285996.229345773*rdxCp2R3[0]*rho[0])*omega*volFac+((12470.76581449591*rdxCp2R4[1]+151831.5737914877*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+187685.0255081635*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-12470.76581449591*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-97895.51164379291*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-871845.0944978701*rdxCp2Sq[0]*rdxCp2Sq[1]-285996.229345773*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-63000.0*rdxCp2[0]*rdxCp2R3[1])-290400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-154800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-10800.0*rdxCp2[0]*rdxCp2R3[1])-66420.0*rdxCp2Sq[0]*rdxCp2Sq[1]-106560.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(12960.0*phiLy[1]-12960.0*phiC[1])*rdxCp2R4[1]+(157788.0*rdxCp2[0]*phiLy[1]-75600.0*rdxCp2[0]*phiLx[1]-535788.0*rdxCp2[0]*phiC[1]+261886.0821044141*rdxCp2[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(465750.0*rdxCp2Sq[0]*phiLy[1]-727155.0*rdxCp2Sq[0]*phiLx[1]-3688425.0*rdxCp2Sq[0]*phiC[1]+1922680.319449907*rdxCp2Sq[0]*bcVals[1]+((-301792.5327108011*phiLy[0])-659547.6270141526*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(195048.0*rdxCp2R3[0]*phiLy[1]-1862820.0*rdxCp2R3[0]*phiLx[1]-7679628.0*rdxCp2R3[0]*phiC[1]+3812313.1094914*rdxCp2R3[0]*bcVals[1]+((-160872.8790069972*phiLy[0])-1745283.675738703*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-681120.0*rdxCp2R4[0]*phiLx[1]-2662560.0*rdxCp2R4[0]*phiC[1]+1286983.032055978*rdxCp2R4[0]*bcVals[1]-643491.516027989*phiLx[0]*rdxCp2R4[0])*omega+12960.0*phiC[1]*rdxCp2R4[1]+535788.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+2662560.0*rdxCp2R4[0]*phiC[1])/(12960.0*rdxCp2R4[1]+535788.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+7679628.0*rdxCp2R3[0]*rdxCp2[1]+2662560.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((144785.5911062975*rdxCp2[0]*rdxCp2Sq[1]+1629679.676638325*rdxCp2Sq[0]*rdxCp2[1]+3860949.096167934*rdxCp2R3[0])*rho[3]+((-25920.0*rdxCp2R3[1])-1006560.0*rdxCp2[0]*rdxCp2Sq[1]-5652000.0*rdxCp2Sq[0]*rdxCp2[1]-8256000.0*rdxCp2R3[0])*rho[2]+(597780.0*rdxCp2[0]*rdxCp2Sq[1]+2352240.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-37412.29744348773*rho[0]*rdxCp2R3[1]-1278253.495985831*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4493632.615156694*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((210444.1731196184*rdxCp2[0]*rdxCp2R3[1]-3016366.481381198*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(16835.53384956948*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]-403117.5049535803*rdxCp2R3[0]*rdxCp2[1]-4289943.440186594*rdxCp2R4[0])*phiLx[3]+((-74824.59488697546*rdxCp2R4[1])-2731097.713374604*rdxCp2[0]*rdxCp2R3[1]-1.26515919188061e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-1.191650955607388e+7*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(94500.0*rdxCp2[0]*rdxCp2R3[1]+1983375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6450000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(9720.0*rdxCp2[0]*rdxCp2R3[1]+63990.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1408860.0*rdxCp2R3[0]*rdxCp2[1]-5263200.0*rdxCp2R4[0])*phiLx[2]+(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0])*phiC[2]+(218700.0*rdxCp2[0]*phiLy[1]+24300.0*rdxCp2[0]*phiLx[1]-224473.7846609264*rdxCp2[0]*bcVals[1]+(98207.28078915528*phiLy[0]+14029.6115413079*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-664200.0*rdxCp2Sq[0]*phiLx[1])-2492594.31717237*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLy[0]-814886.6036909672*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3134700.0*rdxCp2R3[0]*phiLy[1])-2961900.0*rdxCp2R3[0]*phiLx[1]-6590799.732961089*rdxCp2R3[0]*bcVals[1]+(6703036.625291553*phiLy[0]-3407636.758811008*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-64800.0*rdxCp2R4[1])-2678940.0*rdxCp2[0]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]-3.839814e+7*rdxCp2R3[0]*rdxCp2[1]-1.33128e+7*rdxCp2R4[0])*phiC[2]))/(64800.0*rdxCp2R4[1]+2678940.0*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+3.839814e+7*rdxCp2R3[0]*rdxCp2[1]+1.33128e+7*rdxCp2R4[0]); 
  phiC[3] = (((8640.0*rdxCp2R3[1]+253992.0*rdxCp2[0]*rdxCp2Sq[1]+966336.0*rdxCp2Sq[0]*rdxCp2[1]+577920.0*rdxCp2R3[0])*rho[3]+((-17874.76433411081*rdxCp2[0]*rdxCp2Sq[1])-201195.0218072008*rdxCp2Sq[0]*rdxCp2[1]-476660.382242955*rdxCp2R3[0])*rho[2]+(12470.76581449591*rdxCp2R3[1]+89477.74471900817*rdxCp2[0]*rdxCp2Sq[1]+173343.6448214932*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-73800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-290400.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-150000.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-451500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-21600.0*rdxCp2[0]*rdxCp2R3[1])-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2285160.0*rdxCp2R3[0]*rdxCp2[1]-1135200.0*rdxCp2R4[0])*phiLx[3]+((-21600.0*rdxCp2R4[1])-892980.0*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-1.279938e+7*rdxCp2R3[0]*rdxCp2[1]-4437600.0*rdxCp2R4[0])*phiC[3]+((-99600.0*rdxCp2[0]*rdxCp2R3[1])-580800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-688000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(372390.9236273086*rdxCp2R3[0]*rdxCp2[1]-25980.76211353316*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-18706.14872174387*rdxCp2[0]*rdxCp2R3[1])-543205.7742697513*rdxCp2Sq[0]*rdxCp2Sq[1]-2016730.678300898*rdxCp2R3[0]*rdxCp2[1]-1072485.860046648*rdxCp2R4[0])*phiLx[2]+((-155884.5726811989*rdxCp2[0]*phiLy[1])-31176.91453623978*rdxCp2[0]*phiLx[1]+108000.0*rdxCp2[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-687061.2540923841*rdxCp2Sq[0]*phiLy[1])-175759.8556980518*rdxCp2Sq[0]*phiLx[1]+332100.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-469212.5637704087*rdxCp2R3[0]*phiLy[1])-244738.7791094823*rdxCp2R3[0]*phiLx[1]-241200.0*rdxCp2R3[0]*bcVals[1]+(387000.0*phiLy[0]-266400.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0])*phiC[3])/(21600.0*rdxCp2R4[1]+892980.0*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1.279938e+7*rdxCp2R3[0]*rdxCp2[1]+4437600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = -(1.0*(((313632.0*rdxCp2[0]*rdxCp2Sq[1]+79704.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(772189.8192335867*rdxCp2R3[1]+1901292.956078046*rdxCp2[0]*rdxCp2Sq[1]+429306.1131640216*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-599151.015354226*rdxCp2[0]*rdxCp2Sq[1])-170433.7994647775*rdxCp2Sq[0]*rdxCp2[1]-4988.306325798365*rdxCp2R3[0])*rho[1]-1651200.0*rho[0]*rdxCp2R3[1]-4139904.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1057392.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-30240.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-394920.0*rdxCp2[0]*rdxCp2R3[1])-88560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3240.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(417960.0*rdxCp2[0]*rdxCp2R3[1]+784080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+170100.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3219840.0*rdxCp2R4[1])-8373744.0*rdxCp2[0]*rdxCp2R3[1]-2626452.0*rdxCp2Sq[0]*rdxCp2Sq[1]-181440.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-857988.6880373188*rdxCp2R4[1])-2064285.865273508*rdxCp2[0]*rdxCp2R3[1]-386469.0325912283*rdxCp2Sq[0]*rdxCp2Sq[1]+19641.45615783106*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(434356.7733188925*rdxCp2[0]*rdxCp2R3[1]+814839.8383191626*rdxCp2Sq[0]*rdxCp2Sq[1]+176773.1054204795*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(2662560.0*phiC[0]-1052640.0*phiLy[0])*rdxCp2R4[1]+((-454351.5678414679*rdxCp2[0]*phiLy[1])-893738.2167055405*rdxCp2[0]*phiLx[1]-1651200.0*rdxCp2[0]*bcVals[1]+((-2563956.0*phiLy[0])-928800.0*phiLx[0]+7679628.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-108651.5471587956*rdxCp2Sq[0]*phiLy[1])-1772702.040022519*rdxCp2Sq[0]*phiLx[1]-5004704.0*rdxCp2Sq[0]*bcVals[1]+((-532953.0*phiLy[0])-1842246.0*phiLx[0]+3688425.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(1870.614872174387*rdxCp2R3[0]*phiLy[1]-439178.8027671645*rdxCp2R3[0]*phiLx[1]-1303392.0*rdxCp2R3[0]*bcVals[1]+(11340.0*phiLy[0]-456408.0*phiLx[0]+535788.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-12470.76581449591*rdxCp2R4[0]*phiLx[1]-37440.0*rdxCp2R4[0]*bcVals[1]+(12960.0*phiC[0]-12960.0*phiLx[0])*rdxCp2R4[0])*omega-2662560.0*phiC[0]*rdxCp2R4[1]-7679628.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-3688425.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-12960.0*phiC[0]*rdxCp2R4[0]))/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((3860949.096167934*rdxCp2R3[1]+1629679.676638325*rdxCp2[0]*rdxCp2Sq[1]+144785.5911062975*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2352240.0*rdxCp2[0]*rdxCp2Sq[1]+597780.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-8256000.0*rdxCp2R3[1])-5652000.0*rdxCp2[0]*rdxCp2Sq[1]-1006560.0*rdxCp2Sq[0]*rdxCp2[1]-25920.0*rdxCp2R3[0])*rho[1]-4493632.615156694*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1278253.495985831*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-37412.29744348773*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4289943.440186594*rdxCp2R4[1])-403117.5049535803*rdxCp2[0]*rdxCp2R3[1]+291815.9200592043*rdxCp2Sq[0]*rdxCp2Sq[1]+16835.53384956948*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(210444.1731196184*rdxCp2R3[0]*rdxCp2[1]-3016366.481381198*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-6590799.732961089*rdxCp2[0]*rdxCp2R3[1])-2492594.31717237*rdxCp2Sq[0]*rdxCp2Sq[1]-224473.7846609264*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-2961900.0*rdxCp2[0]*rdxCp2R3[1])-664200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(218700.0*rdxCp2R3[0]*rdxCp2[1]-3134700.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(1.33128e+7*phiC[1]-5263200.0*phiLy[1])*rdxCp2R4[1]+((-1408860.0*rdxCp2[0]*phiLy[1])+6450000.0*rdxCp2[0]*phiLx[1]+3.839814e+7*rdxCp2[0]*phiC[1]-1.191650955607388e+7*rdxCp2[0]*bcVals[1]+(6703036.625291553*phiLx[0]-3407636.758811008*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+(63990.0*rdxCp2Sq[0]*phiLy[1]+1983375.0*rdxCp2Sq[0]*phiLx[1]+1.8442125e+7*rdxCp2Sq[0]*phiC[1]-1.26515919188061e+7*rdxCp2Sq[0]*bcVals[1]+(2061183.762277152*phiLx[0]-814886.6036909672*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(9720.0*rdxCp2R3[0]*phiLy[1]+94500.0*rdxCp2R3[0]*phiLx[1]+2678940.0*rdxCp2R3[0]*phiC[1]-2731097.713374604*rdxCp2R3[0]*bcVals[1]+(14029.6115413079*phiLy[0]+98207.28078915528*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+64800.0*rdxCp2R4[0]*phiC[1]-74824.59488697546*rdxCp2R4[0]*bcVals[1])*omega-1.33128e+7*phiC[1]*rdxCp2R4[1]-3.839814e+7*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-1.8442125e+7*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-2678940.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-64800.0*rdxCp2R4[0]*phiC[1]))/(1.33128e+7*rdxCp2R4[1]+3.839814e+7*rdxCp2[0]*rdxCp2R3[1]+1.8442125e+7*rdxCp2Sq[0]*rdxCp2Sq[1]+2678940.0*rdxCp2R3[0]*rdxCp2[1]+64800.0*rdxCp2R4[0]); 
  phiC[2] = (((69337.45792859727*rdxCp2[0]*rdxCp2Sq[1]+35791.09788760327*rdxCp2Sq[0]*rdxCp2[1]+4988.306325798365*rdxCp2R3[0])*rho[3]+(346752.0*rdxCp2R3[1]+928080.0*rdxCp2[0]*rdxCp2Sq[1]+332172.0*rdxCp2Sq[0]*rdxCp2[1]+30240.0*rdxCp2R3[0])*rho[2]+((-116160.0*rdxCp2[0]*rdxCp2Sq[1])-29520.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-285996.229345773*rho[0]*rdxCp2R3[1]-704182.5763252026*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-159002.2641348229*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-97895.51164379291*rdxCp2[0]*rdxCp2R3[1])-70303.9422792207*rdxCp2Sq[0]*rdxCp2Sq[1]-12470.76581449591*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(187685.0255081635*rdxCp2[0]*rdxCp2R3[1]+448168.1464584469*rdxCp2Sq[0]*rdxCp2Sq[1]+151831.5737914877*rdxCp2R3[0]*rdxCp2[1]+12470.76581449591*rdxCp2R4[0])*phiLx[3]+(1286983.032055978*rdxCp2R4[1]+3812313.1094914*rdxCp2[0]*rdxCp2R3[1]+1922680.319449907*rdxCp2Sq[0]*rdxCp2Sq[1]+261886.0821044141*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-681120.0*rdxCp2R4[1])-1862820.0*rdxCp2[0]*rdxCp2R3[1]-727155.0*rdxCp2Sq[0]*rdxCp2Sq[1]-75600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(195048.0*rdxCp2[0]*rdxCp2R3[1]+465750.0*rdxCp2Sq[0]*rdxCp2Sq[1]+157788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiLx[2]+((-2662560.0*rdxCp2R4[1])-7679628.0*rdxCp2[0]*rdxCp2R3[1]-3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]-535788.0*rdxCp2R3[0]*rdxCp2[1]-12960.0*rdxCp2R4[0])*phiC[2]-643491.516027989*phiLy[0]*rdxCp2R4[1]+((-106560.0*rdxCp2[0]*phiLy[1])-154800.0*rdxCp2[0]*phiLx[1]-285996.229345773*rdxCp2[0]*bcVals[1]+((-1745283.675738703*phiLy[0])-160872.8790069972*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-66420.0*rdxCp2Sq[0]*phiLy[1])-290400.0*rdxCp2Sq[0]*phiLx[1]-871845.0944978701*rdxCp2Sq[0]*bcVals[1]+((-659547.6270141526*phiLy[0])-301792.5327108011*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-10800.0*rdxCp2R3[0]*phiLy[1])-63000.0*rdxCp2R3[0]*phiLx[1]-201610.7140010173*rdxCp2R3[0]*bcVals[1]+((-65471.52052610354*phiLy[0])-65471.52052610354*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0])*phiC[2])/(2662560.0*rdxCp2R4[1]+7679628.0*rdxCp2[0]*rdxCp2R3[1]+3688425.0*rdxCp2Sq[0]*rdxCp2Sq[1]+535788.0*rdxCp2R3[0]*rdxCp2[1]+12960.0*rdxCp2R4[0]); 
  phiC[3] = (((577920.0*rdxCp2R3[1]+966336.0*rdxCp2[0]*rdxCp2Sq[1]+253992.0*rdxCp2Sq[0]*rdxCp2[1]+8640.0*rdxCp2R3[0])*rho[3]+(173343.6448214932*rdxCp2[0]*rdxCp2Sq[1]+89477.74471900817*rdxCp2Sq[0]*rdxCp2[1]+12470.76581449591*rdxCp2R3[0])*rho[2]+((-476660.382242955*rdxCp2R3[1])-201195.0218072008*rdxCp2[0]*rdxCp2Sq[1]-17874.76433411081*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-290400.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-73800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1135200.0*rdxCp2R4[1])-2285160.0*rdxCp2[0]*rdxCp2R3[1]-623370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-451500.0*rdxCp2[0]*rdxCp2R3[1])-661125.0*rdxCp2Sq[0]*rdxCp2Sq[1]-150000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-4437600.0*rdxCp2R4[1])-1.279938e+7*rdxCp2[0]*rdxCp2R3[1]-6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892980.0*rdxCp2R3[0]*rdxCp2[1]-21600.0*rdxCp2R4[0])*phiC[3]+((-241200.0*rdxCp2[0]*rdxCp2R3[1])+332100.0*rdxCp2Sq[0]*rdxCp2Sq[1]+108000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-244738.7791094823*rdxCp2[0]*rdxCp2R3[1])-175759.8556980518*rdxCp2Sq[0]*rdxCp2Sq[1]-31176.91453623978*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-469212.5637704087*rdxCp2[0]*rdxCp2R3[1])-687061.2540923841*rdxCp2Sq[0]*rdxCp2Sq[1]-155884.5726811989*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1072485.860046648*phiLy[1]*rdxCp2R4[1]+((-2016730.678300898*rdxCp2[0]*phiLy[1])+372390.9236273086*rdxCp2[0]*phiLx[1]-688000.0*rdxCp2[0]*bcVals[1]+(387000.0*phiLx[0]-266400.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-543205.7742697513*rdxCp2Sq[0]*phiLy[1])-580800.0*rdxCp2Sq[0]*bcVals[1]-166050.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-18706.14872174387*rdxCp2R3[0]*phiLy[1])-25980.76211353316*rdxCp2R3[0]*phiLx[1]-99600.0*rdxCp2R3[0]*bcVals[1]+((-27000.0*phiLy[0])-27000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0])*phiC[3])/(4437600.0*rdxCp2R4[1]+1.279938e+7*rdxCp2[0]*rdxCp2R3[1]+6147375.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892980.0*rdxCp2R3[0]*rdxCp2[1]+21600.0*rdxCp2R4[0]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 
  double *phiUy = phiPrev[3]; 
  double *phiLy = phiPrev[4]; 

  phiC[0] = ((820.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(554.2562584220407*rdxCp2Sq[1]+4416.729559300637*rdxCp2[0]*rdxCp2[1])*rho[2]+(4416.729559300637*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[1]+3360.0*rho[0]*rdxCp2Sq[1]+27351.0*rdxCp2[0]*rho[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((1750.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+1750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(4160.0*rdxCp2R3[1]+33726.0*rdxCp2[0]*rdxCp2Sq[1]+3360.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1385.640646055102*rdxCp2R3[1]+11353.59304361399*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+1818.653347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(1440.0*phiLy[0]-1440.0*phiC[0])*rdxCp2R3[1]+(1818.653347947321*rdxCp2[0]*phiLy[1]+1818.653347947321*rdxCp2[0]*phiLx[1]+3360.0*rdxCp2[0]*bcVals[1]+(11799.0*phiLy[0]+1890.0*phiLx[0]-13689.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2Sq[0]*phiLy[1]+11353.59304361399*rdxCp2Sq[0]*phiLx[1]+33726.0*rdxCp2Sq[0]*bcVals[1]+(1890.0*phiLy[0]+11799.0*phiLx[0]-13689.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+1385.640646055102*rdxCp2R3[0]*phiLx[1]+4160.0*rdxCp2R3[0]*bcVals[1]+(1440.0*phiLx[0]-1440.0*phiC[0])*rdxCp2R3[0])*omega+1440.0*phiC[0]*rdxCp2R3[1]+13689.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+13689.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+1440.0*phiC[0]*rdxCp2R3[0])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[1] = (((554.2562584220407*rdxCp2Sq[1]+297.9127389018469*rdxCp2[0]*rdxCp2[1])*rho[3]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(3360.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+576.0*rdxCp2Sq[0])*rho[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1]+831.384387633061*rdxCp2Sq[0]*rho[0])*omega*volFac+((1385.640646055102*rdxCp2R3[1]+2563.435195201938*rdxCp2[0]*rdxCp2Sq[1]+311.7691453623978*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(8400.446416709054*rdxCp2[0]*rdxCp2Sq[1]+831.384387633061*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2625.0*rdxCp2[0]*rdxCp2Sq[1]+450.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(450.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(1440.0*phiLy[1]-1440.0*phiC[1])*rdxCp2R3[1]+(2664.0*rdxCp2[0]*phiLy[1]-2625.0*rdxCp2[0]*phiLx[1]-13689.0*rdxCp2[0]*phiC[1]+4849.742261192856*rdxCp2[0]*bcVals[1]+(2727.980021920981*phiLy[0]-2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(324.0*rdxCp2Sq[0]*phiLy[1]-450.0*rdxCp2Sq[0]*phiLx[1]-13689.0*rdxCp2Sq[0]*phiC[1]+14081.57306553497*rdxCp2Sq[0]*bcVals[1]+(467.6537180435967*phiLy[0]-467.6537180435967*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-1440.0*rdxCp2R3[0]*phiC[1]+1662.768775266122*rdxCp2R3[0]*bcVals[1])*omega+1440.0*phiC[1]*rdxCp2R3[1]+13689.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+1440.0*rdxCp2R3[0]*phiC[1])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[2] = (((297.9127389018469*rdxCp2[0]*rdxCp2[1]+554.2562584220407*rdxCp2Sq[0])*rho[3]+(576.0*rdxCp2Sq[1]+5166.0*rdxCp2[0]*rdxCp2[1]+3360.0*rdxCp2Sq[0])*rho[2]+1230.0*rdxCp2[0]*rdxCp2[1]*rho[1]+831.384387633061*rho[0]*rdxCp2Sq[1]+6625.094338950954*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(311.7691453623978*rdxCp2[0]*rdxCp2Sq[1]+2563.435195201938*rdxCp2Sq[0]*rdxCp2[1]+1385.640646055102*rdxCp2R3[0])*phiLx[3]+(1662.768775266122*rdxCp2R3[1]+14081.57306553497*rdxCp2[0]*rdxCp2Sq[1]+4849.742261192856*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-450.0*rdxCp2[0]*rdxCp2Sq[1])-2625.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(324.0*rdxCp2[0]*rdxCp2Sq[1]+2664.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiLx[2]+((-1440.0*rdxCp2R3[1])-13689.0*rdxCp2[0]*rdxCp2Sq[1]-13689.0*rdxCp2Sq[0]*rdxCp2[1]-1440.0*rdxCp2R3[0])*phiC[2]+(450.0*rdxCp2[0]*phiLy[1]+450.0*rdxCp2[0]*phiLx[1]+831.384387633061*rdxCp2[0]*bcVals[1]+(467.6537180435967*phiLx[0]-467.6537180435967*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-450.0*rdxCp2Sq[0]*phiLy[1])+2625.0*rdxCp2Sq[0]*phiLx[1]+8400.446416709054*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiLx[0]-2727.980021920981*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0])*phiC[2])/(1440.0*rdxCp2R3[1]+13689.0*rdxCp2[0]*rdxCp2Sq[1]+13689.0*rdxCp2Sq[0]*rdxCp2[1]+1440.0*rdxCp2R3[0]); 
  phiC[3] = (((960.0*rdxCp2Sq[1]+6116.0*rdxCp2[0]*rdxCp2[1]+960.0*rdxCp2Sq[0])*rho[3]+(744.7818472546172*rdxCp2[0]*rdxCp2[1]+1385.640646055102*rdxCp2Sq[0])*rho[2]+(1385.640646055102*rdxCp2Sq[1]+744.7818472546172*rdxCp2[0]*rdxCp2[1])*rho[1]+3075.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4375.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-4375.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-2400.0*rdxCp2R3[1])-22815.0*rdxCp2[0]*rdxCp2Sq[1]-22815.0*rdxCp2Sq[0]*rdxCp2[1]-2400.0*rdxCp2R3[0])*phiC[3]+(4150.0*rdxCp2[0]*rdxCp2Sq[1]+2000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1082.531754730548*rdxCp2[0]*rdxCp2Sq[1]-1082.531754730548*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-779.4228634059946*rdxCp2[0]*rdxCp2Sq[1])-4546.633369868302*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-4546.633369868302*rdxCp2[0]*phiLy[1])-1082.531754730548*rdxCp2[0]*phiLx[1]+2000.0*rdxCp2[0]*bcVals[1]+(1125.0*phiLy[0]-1125.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-779.4228634059946*rdxCp2Sq[0]*phiLy[1])+1082.531754730548*rdxCp2Sq[0]*phiLx[1]+4150.0*rdxCp2Sq[0]*bcVals[1]+(1125.0*phiLx[0]-1125.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0])*phiC[3])/(2400.0*rdxCp2R3[1]+22815.0*rdxCp2[0]*rdxCp2Sq[1]+22815.0*rdxCp2Sq[0]*rdxCp2[1]+2400.0*rdxCp2R3[0]); 

}

void MGpoissonResidue2xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 
  double *dxUy = dx[3]; 
  double *dxLy = dx[4]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdxCp2[2]; 
  double rdxLx[2]; 
  double rdxUx[2]; 
  double rdxLxSq[2]; 
  double rdxUxSq[2]; 
  double rdxLxCu[2]; 
  double rdxUxCu[2]; 
  double rdxLxR4[2]; 
  double rdxUxR4[2]; 
  double rdxLy[2]; 
  double rdxUy[2]; 
  double rdxLySq[2]; 
  double rdxUySq[2]; 
  double rdxLyCu[2]; 
  double rdxUyCu[2]; 
  double rdxLyR4[2]; 
  double rdxUyR4[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = volFac*4.0/(dxLx[1]*dxLx[1]); 
  rdxUx[1]   = volFac*4.0/(dxUx[1]*dxUx[1]); 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = volFac*4.0/(dxLy[0]*dxLy[0]); 
  rdxUy[0]   = volFac*4.0/(dxUy[0]*dxUy[0]); 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = volFac*4.0/(dxLy[1]*dxLy[1]); 
  rdxUy[1]   = volFac*4.0/(dxUy[1]*dxUy[1]); 
  rdxLySq[1] = rdxLy[1]*rdxLy[1]; 
  rdxUySq[1] = rdxUy[1]*rdxUy[1]; 
  rdxLyCu[1] = rdxLySq[1]*rdxLy[1]; 
  rdxUyCu[1] = rdxUySq[1]*rdxUy[1]; 
  rdxLyR4[1] = rdxLyCu[1]*rdxLy[1]; 
  rdxUyR4[1] = rdxUyCu[1]*rdxUy[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxUy[1]*phiUy[2]+8.660254037844386*rdxLy[1]*phiLy[2]+(8.660254037844386*rdxLy[1]-8.660254037844386*rdxUy[1])*phiC[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdxUy[1]+(9.0*phiLy[0]-9.0*phiC[0])*rdxLy[1]-8.660254037844386*rdxUx[0]*phiUx[1]+8.660254037844386*rdxLx[0]*phiLx[1]+(8.660254037844386*rdxLx[0]-8.660254037844386*rdxUx[0])*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0])*rdxUx[0]+(9.0*phiLx[0]-9.0*phiC[0])*rdxLx[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdxUy[1]*phiUy[3]+8.660254037844386*rdxLy[1]*phiLy[3]+(8.660254037844386*rdxLy[1]-8.660254037844386*rdxUy[1])*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdxUy[1]+(9.0*phiLy[1]-9.0*phiC[1])*rdxLy[1]-7.0*rdxUx[0]*phiUx[1]-7.0*rdxLx[0]*phiLx[1]+((-23.0*rdxUx[0])-23.0*rdxLx[0])*phiC[1]+(8.660254037844386*phiUx[0]-22.5166604983954*phiC[0])*rdxUx[0]+(22.5166604983954*phiC[0]-8.660254037844386*phiLx[0])*rdxLx[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdxUx[0]*phiUx[3]+8.660254037844386*rdxLx[0]*phiLx[3]+(8.660254037844386*rdxLx[0]-8.660254037844386*rdxUx[0])*phiC[3]-7.0*rdxUy[1]*phiUy[2]+9.0*rdxUx[0]*phiUx[2]-7.0*rdxLy[1]*phiLy[2]+9.0*rdxLx[0]*phiLx[2]+((-23.0*rdxUy[1])-23.0*rdxLy[1]-9.0*rdxUx[0]-9.0*rdxLx[0])*phiC[2]+(8.660254037844386*phiUy[0]-22.5166604983954*phiC[0])*rdxUy[1]+(22.5166604983954*phiC[0]-8.660254037844386*phiLy[0])*rdxLy[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdxUy[1]*phiUy[3]-7.0*rdxUx[0]*phiUx[3]-7.0*rdxLy[1]*phiLy[3]-7.0*rdxLx[0]*phiLx[3]+((-23.0*rdxUy[1])-23.0*rdxLy[1]-23.0*rdxUx[0]-23.0*rdxLx[0])*phiC[3]+8.660254037844386*rdxUx[0]*phiUx[2]-8.660254037844386*rdxLx[0]*phiLx[2]+(22.5166604983954*rdxLx[0]-22.5166604983954*rdxUx[0])*phiC[2]+(8.660254037844386*phiUy[1]-22.5166604983954*phiC[1])*rdxUy[1]+(22.5166604983954*phiC[1]-8.660254037844386*phiLy[1])*rdxLy[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]+8.660254037844386*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiUx[1]+93.53074360871933*rdxCp2[0]*phiC[1]+((-6.0*phiUx[0])-42.0*phiC[0]+96.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdxCp2[1]*phiUy[3]+8.660254037844386*rdxCp2[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1]-18.0*phiC[1])*rdxCp2[1]-40.0*rdxCp2[0]*phiUx[1]-200.0*rdxCp2[0]*phiC[1]+(34.64101615137754*phiUx[0]+34.64101615137754*phiC[0]-138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+10.39230484541326*rdxCp2[0]*phiUx[3]+93.53074360871933*rdxCp2[0]*phiC[3]-7.0*rdxCp2[1]*phiUy[2]-6.0*rdxCp2[0]*phiUx[2]-7.0*rdxCp2[1]*phiLy[2]+((-46.0*rdxCp2[1])-42.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdxCp2[1]*phiUy[3]-40.0*rdxCp2[0]*phiUx[3]-7.0*rdxCp2[1]*phiLy[3]+((-46.0*rdxCp2[1])-200.0*rdxCp2[0])*phiC[3]+34.64101615137754*rdxCp2[0]*phiUx[2]+34.64101615137754*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]+8.660254037844386*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]-8.660254037844386*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac-51.96152422706631*rdxCp2[1]*phiUy[3]+51.96152422706631*rdxCp2[1]*phiLy[3]+(54.0*phiUy[1]+54.0*phiLy[1]-108.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiUx[1]-315.0*rdxCp2[0]*phiC[1]+(77.94228634059945*phiUx[0]-77.94228634059945*phiC[0]+138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdxCp2[0]*phiUx[3]-8.660254037844386*rdxCp2[0]*phiC[3]-7.0*rdxCp2[1]*phiUy[2]+9.0*rdxCp2[0]*phiUx[2]-7.0*rdxCp2[1]*phiLy[2]+((-46.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-14.0*rdxCp2[1]*phiUy[3]-25.0*rdxCp2[0]*phiUx[3]-14.0*rdxCp2[1]*phiLy[3]+((-92.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]+8.660254037844386*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-93.53074360871933*rdxCp2[0]*phiC[1]+96.0*rdxCp2[0]*bcVals[1]+((-6.0*phiLx[0])-42.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdxCp2[1]*phiUy[3]+8.660254037844386*rdxCp2[1]*phiLy[3]+(9.0*phiUy[1]+9.0*phiLy[1]-18.0*phiC[1])*rdxCp2[1]-40.0*rdxCp2[0]*phiLx[1]-200.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+((-34.64101615137754*phiLx[0])-34.64101615137754*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-10.39230484541326*rdxCp2[0]*phiLx[3]-93.53074360871933*rdxCp2[0]*phiC[3]-7.0*rdxCp2[1]*phiUy[2]-7.0*rdxCp2[1]*phiLy[2]-6.0*rdxCp2[0]*phiLx[2]+((-46.0*rdxCp2[1])-42.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdxCp2[1]*phiUy[3]-7.0*rdxCp2[1]*phiLy[3]-40.0*rdxCp2[0]*phiLx[3]+((-46.0*rdxCp2[1])-200.0*rdxCp2[0])*phiC[3]-34.64101615137754*rdxCp2[0]*phiLx[2]-34.64101615137754*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]+8.660254037844386*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+8.660254037844386*rdxCp2[0]*phiC[1]+16.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac-51.96152422706631*rdxCp2[1]*phiUy[3]+51.96152422706631*rdxCp2[1]*phiLy[3]+(54.0*phiUy[1]+54.0*phiLy[1]-108.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiLx[1]-315.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiC[0]-77.94228634059945*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+8.660254037844386*rdxCp2[0]*phiLx[3]+8.660254037844386*rdxCp2[0]*phiC[3]-7.0*rdxCp2[1]*phiUy[2]-7.0*rdxCp2[1]*phiLy[2]+9.0*rdxCp2[0]*phiLx[2]+((-46.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-14.0*rdxCp2[1]*phiUy[3]-14.0*rdxCp2[1]*phiLy[3]-25.0*rdxCp2[0]*phiLx[3]+((-92.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+10.39230484541326*rdxCp2[1]*phiUy[2]+93.53074360871933*rdxCp2[1]*phiC[2]+96.0*rdxCp2[1]*bcVals[2]+((-6.0*phiUy[0])-42.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+10.39230484541326*rdxCp2[1]*phiUy[3]+93.53074360871933*rdxCp2[1]*phiC[3]+((-6.0*phiUy[1])-42.0*phiC[1])*rdxCp2[1]-7.0*rdxCp2[0]*phiUx[1]-7.0*rdxCp2[0]*phiLx[1]-46.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdxCp2[0]*phiUx[3]+8.660254037844386*rdxCp2[0]*phiLx[3]-40.0*rdxCp2[1]*phiUy[2]+9.0*rdxCp2[0]*phiUx[2]+9.0*rdxCp2[0]*phiLx[2]+((-200.0*rdxCp2[1])-18.0*rdxCp2[0])*phiC[2]-138.5640646055102*rdxCp2[1]*bcVals[2]+(34.64101615137754*phiUy[0]+34.64101615137754*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-40.0*rdxCp2[1]*phiUy[3]-7.0*rdxCp2[0]*phiUx[3]-7.0*rdxCp2[0]*phiLx[3]+((-200.0*rdxCp2[1])-46.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]-8.660254037844386*rdxCp2[0]*phiLx[2]+(34.64101615137754*phiUy[1]+34.64101615137754*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]-8.660254037844386*rdxCp2[1]*phiC[2]-16.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdxCp2[1]*phiUy[3]-8.660254037844386*rdxCp2[1]*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-7.0*rdxCp2[0]*phiUx[1]-7.0*rdxCp2[0]*phiLx[1]-46.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac-51.96152422706631*rdxCp2[0]*phiUx[3]+51.96152422706631*rdxCp2[0]*phiLx[3]-75.0*rdxCp2[1]*phiUy[2]+54.0*rdxCp2[0]*phiUx[2]+54.0*rdxCp2[0]*phiLx[2]+((-315.0*rdxCp2[1])-108.0*rdxCp2[0])*phiC[2]+138.5640646055102*rdxCp2[1]*bcVals[2]+(77.94228634059945*phiUy[0]-77.94228634059945*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiUy[3]-14.0*rdxCp2[0]*phiUx[3]-14.0*rdxCp2[0]*phiLx[3]+((-105.0*rdxCp2[1])-92.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiLx[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+96.0*rdxCp2[1]*bcVals[3]-10.39230484541326*rdxCp2[1]*phiLy[2]-93.53074360871933*rdxCp2[1]*phiC[2]+((-6.0*phiLy[0])-42.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-10.39230484541326*rdxCp2[1]*phiLy[3]-93.53074360871933*rdxCp2[1]*phiC[3]+((-6.0*phiLy[1])-42.0*phiC[1])*rdxCp2[1]-7.0*rdxCp2[0]*phiUx[1]-7.0*rdxCp2[0]*phiLx[1]-46.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdxCp2[0]*phiUx[3]+8.660254037844386*rdxCp2[0]*phiLx[3]+138.5640646055102*rdxCp2[1]*bcVals[3]+9.0*rdxCp2[0]*phiUx[2]-40.0*rdxCp2[1]*phiLy[2]+9.0*rdxCp2[0]*phiLx[2]+((-200.0*rdxCp2[1])-18.0*rdxCp2[0])*phiC[2]+((-34.64101615137754*phiLy[0])-34.64101615137754*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.0625*(16.0*rho[3]*volFac-7.0*rdxCp2[0]*phiUx[3]-40.0*rdxCp2[1]*phiLy[3]-7.0*rdxCp2[0]*phiLx[3]+((-200.0*rdxCp2[1])-46.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]-8.660254037844386*rdxCp2[0]*phiLx[2]+((-34.64101615137754*phiLy[1])-34.64101615137754*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdxCp2[1]*bcVals[3]+8.660254037844386*rdxCp2[1]*phiLy[2]+8.660254037844386*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+8.660254037844386*rdxCp2[1]*phiLy[3]+8.660254037844386*rdxCp2[1]*phiC[3]+(9.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-7.0*rdxCp2[0]*phiUx[1]-7.0*rdxCp2[0]*phiLx[1]-46.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac-51.96152422706631*rdxCp2[0]*phiUx[3]+51.96152422706631*rdxCp2[0]*phiLx[3]+138.5640646055102*rdxCp2[1]*bcVals[3]+54.0*rdxCp2[0]*phiUx[2]-75.0*rdxCp2[1]*phiLy[2]+54.0*rdxCp2[0]*phiLx[2]+((-315.0*rdxCp2[1])-108.0*rdxCp2[0])*phiC[2]+(77.94228634059945*phiC[0]-77.94228634059945*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-14.0*rdxCp2[0]*phiUx[3]-25.0*rdxCp2[1]*phiLy[3]-14.0*rdxCp2[0]*phiLx[3]+((-105.0*rdxCp2[1])-92.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiLx[2]+(25.98076211353316*phiC[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.125*(8.0*rho[0]*volFac+5.196152422706631*rdxCp2[1]*phiUy[2]+46.76537180435967*rdxCp2[1]*phiC[2]+48.0*rdxCp2[1]*bcVals[2]+((-3.0*phiUy[0])-21.0*phiC[0])*rdxCp2[1]+5.196152422706631*rdxCp2[0]*phiUx[1]+46.76537180435967*rdxCp2[0]*phiC[1]+((-3.0*phiUx[0])-21.0*phiC[0]+48.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.125*(8.0*rho[1]*volFac+5.196152422706631*rdxCp2[1]*phiUy[3]+46.76537180435967*rdxCp2[1]*phiC[3]+((-3.0*phiUy[1])-21.0*phiC[1])*rdxCp2[1]-20.0*rdxCp2[0]*phiUx[1]-100.0*rdxCp2[0]*phiC[1]+(17.32050807568877*phiUx[0]+17.32050807568877*phiC[0]-69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.125*(8.0*rho[2]*volFac+5.196152422706631*rdxCp2[0]*phiUx[3]+46.76537180435967*rdxCp2[0]*phiC[3]-20.0*rdxCp2[1]*phiUy[2]-3.0*rdxCp2[0]*phiUx[2]+((-100.0*rdxCp2[1])-21.0*rdxCp2[0])*phiC[2]-69.28203230275508*rdxCp2[1]*bcVals[2]+(17.32050807568877*phiUy[0]+17.32050807568877*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[0]*phiUx[3]+((-50.0*rdxCp2[1])-50.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]+8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]+8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]-8.660254037844386*rdxCp2[1]*phiC[2]-16.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiUx[1]+93.53074360871933*rdxCp2[0]*phiC[1]+((-6.0*phiUx[0])-42.0*phiC[0]+96.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdxCp2[1]*phiUy[3]-8.660254037844386*rdxCp2[1]*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-40.0*rdxCp2[0]*phiUx[1]-200.0*rdxCp2[0]*phiC[1]+(34.64101615137754*phiUx[0]+34.64101615137754*phiC[0]-138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac+62.35382907247956*rdxCp2[0]*phiUx[3]+561.1844616523159*rdxCp2[0]*phiC[3]-75.0*rdxCp2[1]*phiUy[2]-36.0*rdxCp2[0]*phiUx[2]+((-315.0*rdxCp2[1])-252.0*rdxCp2[0])*phiC[2]+138.5640646055102*rdxCp2[1]*bcVals[2]+(77.94228634059945*phiUy[0]-77.94228634059945*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiUy[3]-80.0*rdxCp2[0]*phiUx[3]+((-105.0*rdxCp2[1])-400.0*rdxCp2[0])*phiC[3]+69.28203230275508*rdxCp2[0]*phiUx[2]+69.28203230275508*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+10.39230484541326*rdxCp2[1]*phiUy[2]+93.53074360871933*rdxCp2[1]*phiC[2]+96.0*rdxCp2[1]*bcVals[2]+((-6.0*phiUy[0])-42.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]-8.660254037844386*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac+62.35382907247956*rdxCp2[1]*phiUy[3]+561.1844616523159*rdxCp2[1]*phiC[3]+((-36.0*phiUy[1])-252.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiUx[1]-315.0*rdxCp2[0]*phiC[1]+(77.94228634059945*phiUx[0]-77.94228634059945*phiC[0]+138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdxCp2[0]*phiUx[3]-8.660254037844386*rdxCp2[0]*phiC[3]-40.0*rdxCp2[1]*phiUy[2]+9.0*rdxCp2[0]*phiUx[2]+((-200.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]-138.5640646055102*rdxCp2[1]*bcVals[2]+(34.64101615137754*phiUy[0]+34.64101615137754*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-80.0*rdxCp2[1]*phiUy[3]-25.0*rdxCp2[0]*phiUx[3]+((-400.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+(69.28203230275508*phiUy[1]+69.28203230275508*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]-8.660254037844386*rdxCp2[1]*phiC[2]-16.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]-8.660254037844386*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac-51.96152422706631*rdxCp2[1]*phiUy[3]-51.96152422706631*rdxCp2[1]*phiC[3]+(54.0*phiUy[1]-54.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiUx[1]-315.0*rdxCp2[0]*phiC[1]+(77.94228634059945*phiUx[0]-77.94228634059945*phiC[0]+138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac-51.96152422706631*rdxCp2[0]*phiUx[3]-51.96152422706631*rdxCp2[0]*phiC[3]-75.0*rdxCp2[1]*phiUy[2]+54.0*rdxCp2[0]*phiUx[2]+((-315.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+138.5640646055102*rdxCp2[1]*bcVals[2]+(77.94228634059945*phiUy[0]-77.94228634059945*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiUy[3]-25.0*rdxCp2[0]*phiUx[3]+((-105.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.125*(8.0*rho[0]*volFac+48.0*rdxCp2[1]*bcVals[3]-5.196152422706631*rdxCp2[1]*phiLy[2]-46.76537180435967*rdxCp2[1]*phiC[2]+((-3.0*phiLy[0])-21.0*phiC[0])*rdxCp2[1]+5.196152422706631*rdxCp2[0]*phiUx[1]+46.76537180435967*rdxCp2[0]*phiC[1]+((-3.0*phiUx[0])-21.0*phiC[0]+48.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.125*(8.0*rho[1]*volFac-5.196152422706631*rdxCp2[1]*phiLy[3]-46.76537180435967*rdxCp2[1]*phiC[3]+((-3.0*phiLy[1])-21.0*phiC[1])*rdxCp2[1]-20.0*rdxCp2[0]*phiUx[1]-100.0*rdxCp2[0]*phiC[1]+(17.32050807568877*phiUx[0]+17.32050807568877*phiC[0]-69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.125*(8.0*rho[2]*volFac+5.196152422706631*rdxCp2[0]*phiUx[3]+46.76537180435967*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]-3.0*rdxCp2[0]*phiUx[2]-20.0*rdxCp2[1]*phiLy[2]+((-100.0*rdxCp2[1])-21.0*rdxCp2[0])*phiC[2]+((-17.32050807568877*phiLy[0])-17.32050807568877*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[0]*phiUx[3]-10.0*rdxCp2[1]*phiLy[3]+((-50.0*rdxCp2[1])-50.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]+8.660254037844386*rdxCp2[0]*phiC[2]+((-8.660254037844386*phiLy[1])-8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdxCp2[1]*bcVals[3]+8.660254037844386*rdxCp2[1]*phiLy[2]+8.660254037844386*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiUx[1]+93.53074360871933*rdxCp2[0]*phiC[1]+((-6.0*phiUx[0])-42.0*phiC[0]+96.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+8.660254037844386*rdxCp2[1]*phiLy[3]+8.660254037844386*rdxCp2[1]*phiC[3]+(9.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-40.0*rdxCp2[0]*phiUx[1]-200.0*rdxCp2[0]*phiC[1]+(34.64101615137754*phiUx[0]+34.64101615137754*phiC[0]-138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac+62.35382907247956*rdxCp2[0]*phiUx[3]+561.1844616523159*rdxCp2[0]*phiC[3]+138.5640646055102*rdxCp2[1]*bcVals[3]-36.0*rdxCp2[0]*phiUx[2]-75.0*rdxCp2[1]*phiLy[2]+((-315.0*rdxCp2[1])-252.0*rdxCp2[0])*phiC[2]+(77.94228634059945*phiC[0]-77.94228634059945*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-80.0*rdxCp2[0]*phiUx[3]-25.0*rdxCp2[1]*phiLy[3]+((-105.0*rdxCp2[1])-400.0*rdxCp2[0])*phiC[3]+69.28203230275508*rdxCp2[0]*phiUx[2]+69.28203230275508*rdxCp2[0]*phiC[2]+(25.98076211353316*phiC[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+96.0*rdxCp2[1]*bcVals[3]-10.39230484541326*rdxCp2[1]*phiLy[2]-93.53074360871933*rdxCp2[1]*phiC[2]+((-6.0*phiLy[0])-42.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]-8.660254037844386*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac-62.35382907247956*rdxCp2[1]*phiLy[3]-561.1844616523159*rdxCp2[1]*phiC[3]+((-36.0*phiLy[1])-252.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiUx[1]-315.0*rdxCp2[0]*phiC[1]+(77.94228634059945*phiUx[0]-77.94228634059945*phiC[0]+138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-8.660254037844386*rdxCp2[0]*phiUx[3]-8.660254037844386*rdxCp2[0]*phiC[3]+138.5640646055102*rdxCp2[1]*bcVals[3]+9.0*rdxCp2[0]*phiUx[2]-40.0*rdxCp2[1]*phiLy[2]+((-200.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+((-34.64101615137754*phiLy[0])-34.64101615137754*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[0]*phiUx[3]-80.0*rdxCp2[1]*phiLy[3]+((-400.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+((-69.28203230275508*phiLy[1])-69.28203230275508*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdxCp2[1]*bcVals[3]+8.660254037844386*rdxCp2[1]*phiLy[2]+8.660254037844386*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-8.660254037844386*rdxCp2[0]*phiUx[1]-8.660254037844386*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0]-16.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac+51.96152422706631*rdxCp2[1]*phiLy[3]+51.96152422706631*rdxCp2[1]*phiC[3]+(54.0*phiLy[1]-54.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiUx[1]-315.0*rdxCp2[0]*phiC[1]+(77.94228634059945*phiUx[0]-77.94228634059945*phiC[0]+138.5640646055102*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac-51.96152422706631*rdxCp2[0]*phiUx[3]-51.96152422706631*rdxCp2[0]*phiC[3]+138.5640646055102*rdxCp2[1]*bcVals[3]+54.0*rdxCp2[0]*phiUx[2]-75.0*rdxCp2[1]*phiLy[2]+((-315.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+(77.94228634059945*phiC[0]-77.94228634059945*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[0]*phiUx[3]-25.0*rdxCp2[1]*phiLy[3]+((-105.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+(25.98076211353316*phiC[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.125*(8.0*rho[0]*volFac+5.196152422706631*rdxCp2[1]*phiUy[2]+46.76537180435967*rdxCp2[1]*phiC[2]+48.0*rdxCp2[1]*bcVals[2]+((-3.0*phiUy[0])-21.0*phiC[0])*rdxCp2[1]-5.196152422706631*rdxCp2[0]*phiLx[1]-46.76537180435967*rdxCp2[0]*phiC[1]+48.0*rdxCp2[0]*bcVals[1]+((-3.0*phiLx[0])-21.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.125*(8.0*rho[1]*volFac+5.196152422706631*rdxCp2[1]*phiUy[3]+46.76537180435967*rdxCp2[1]*phiC[3]+((-3.0*phiUy[1])-21.0*phiC[1])*rdxCp2[1]-20.0*rdxCp2[0]*phiLx[1]-100.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+((-17.32050807568877*phiLx[0])-17.32050807568877*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.125*(8.0*rho[2]*volFac-5.196152422706631*rdxCp2[0]*phiLx[3]-46.76537180435967*rdxCp2[0]*phiC[3]-20.0*rdxCp2[1]*phiUy[2]-3.0*rdxCp2[0]*phiLx[2]+((-100.0*rdxCp2[1])-21.0*rdxCp2[0])*phiC[2]-69.28203230275508*rdxCp2[1]*bcVals[2]+(17.32050807568877*phiUy[0]+17.32050807568877*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[0]*phiLx[3]+((-50.0*rdxCp2[1])-50.0*rdxCp2[0])*phiC[3]-8.660254037844386*rdxCp2[0]*phiLx[2]-8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]+8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]-8.660254037844386*rdxCp2[1]*phiC[2]-16.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-93.53074360871933*rdxCp2[0]*phiC[1]+96.0*rdxCp2[0]*bcVals[1]+((-6.0*phiLx[0])-42.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-8.660254037844386*rdxCp2[1]*phiUy[3]-8.660254037844386*rdxCp2[1]*phiC[3]+(9.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-40.0*rdxCp2[0]*phiLx[1]-200.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+((-34.64101615137754*phiLx[0])-34.64101615137754*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac-62.35382907247956*rdxCp2[0]*phiLx[3]-561.1844616523159*rdxCp2[0]*phiC[3]-75.0*rdxCp2[1]*phiUy[2]-36.0*rdxCp2[0]*phiLx[2]+((-315.0*rdxCp2[1])-252.0*rdxCp2[0])*phiC[2]+138.5640646055102*rdxCp2[1]*bcVals[2]+(77.94228634059945*phiUy[0]-77.94228634059945*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiUy[3]-80.0*rdxCp2[0]*phiLx[3]+((-105.0*rdxCp2[1])-400.0*rdxCp2[0])*phiC[3]-69.28203230275508*rdxCp2[0]*phiLx[2]-69.28203230275508*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+10.39230484541326*rdxCp2[1]*phiUy[2]+93.53074360871933*rdxCp2[1]*phiC[2]+96.0*rdxCp2[1]*bcVals[2]+((-6.0*phiUy[0])-42.0*phiC[0])*rdxCp2[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+8.660254037844386*rdxCp2[0]*phiC[1]+16.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac+62.35382907247956*rdxCp2[1]*phiUy[3]+561.1844616523159*rdxCp2[1]*phiC[3]+((-36.0*phiUy[1])-252.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiLx[1]-315.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiC[0]-77.94228634059945*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+8.660254037844386*rdxCp2[0]*phiLx[3]+8.660254037844386*rdxCp2[0]*phiC[3]-40.0*rdxCp2[1]*phiUy[2]+9.0*rdxCp2[0]*phiLx[2]+((-200.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]-138.5640646055102*rdxCp2[1]*bcVals[2]+(34.64101615137754*phiUy[0]+34.64101615137754*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-80.0*rdxCp2[1]*phiUy[3]-25.0*rdxCp2[0]*phiLx[3]+((-400.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+(69.28203230275508*phiUy[1]+69.28203230275508*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxCp2[1]*phiUy[2]-8.660254037844386*rdxCp2[1]*phiC[2]-16.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+8.660254037844386*rdxCp2[0]*phiC[1]+16.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac-51.96152422706631*rdxCp2[1]*phiUy[3]-51.96152422706631*rdxCp2[1]*phiC[3]+(54.0*phiUy[1]-54.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiLx[1]-315.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiC[0]-77.94228634059945*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac+51.96152422706631*rdxCp2[0]*phiLx[3]+51.96152422706631*rdxCp2[0]*phiC[3]-75.0*rdxCp2[1]*phiUy[2]+54.0*rdxCp2[0]*phiLx[2]+((-315.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+138.5640646055102*rdxCp2[1]*bcVals[2]+(77.94228634059945*phiUy[0]-77.94228634059945*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiUy[3]-25.0*rdxCp2[0]*phiLx[3]+((-105.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.125*(8.0*rho[0]*volFac+48.0*rdxCp2[1]*bcVals[3]-5.196152422706631*rdxCp2[1]*phiLy[2]-46.76537180435967*rdxCp2[1]*phiC[2]+((-3.0*phiLy[0])-21.0*phiC[0])*rdxCp2[1]-5.196152422706631*rdxCp2[0]*phiLx[1]-46.76537180435967*rdxCp2[0]*phiC[1]+48.0*rdxCp2[0]*bcVals[1]+((-3.0*phiLx[0])-21.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.125*(8.0*rho[1]*volFac-5.196152422706631*rdxCp2[1]*phiLy[3]-46.76537180435967*rdxCp2[1]*phiC[3]+((-3.0*phiLy[1])-21.0*phiC[1])*rdxCp2[1]-20.0*rdxCp2[0]*phiLx[1]-100.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+((-17.32050807568877*phiLx[0])-17.32050807568877*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.125*(8.0*rho[2]*volFac-5.196152422706631*rdxCp2[0]*phiLx[3]-46.76537180435967*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]-20.0*rdxCp2[1]*phiLy[2]-3.0*rdxCp2[0]*phiLx[2]+((-100.0*rdxCp2[1])-21.0*rdxCp2[0])*phiC[2]+((-17.32050807568877*phiLy[0])-17.32050807568877*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[1]*phiLy[3]-10.0*rdxCp2[0]*phiLx[3]+((-50.0*rdxCp2[1])-50.0*rdxCp2[0])*phiC[3]-8.660254037844386*rdxCp2[0]*phiLx[2]-8.660254037844386*rdxCp2[0]*phiC[2]+((-8.660254037844386*phiLy[1])-8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdxCp2[1]*bcVals[3]+8.660254037844386*rdxCp2[1]*phiLy[2]+8.660254037844386*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-93.53074360871933*rdxCp2[0]*phiC[1]+96.0*rdxCp2[0]*bcVals[1]+((-6.0*phiLx[0])-42.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac+8.660254037844386*rdxCp2[1]*phiLy[3]+8.660254037844386*rdxCp2[1]*phiC[3]+(9.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-40.0*rdxCp2[0]*phiLx[1]-200.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+((-34.64101615137754*phiLx[0])-34.64101615137754*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac-62.35382907247956*rdxCp2[0]*phiLx[3]-561.1844616523159*rdxCp2[0]*phiC[3]+138.5640646055102*rdxCp2[1]*bcVals[3]-75.0*rdxCp2[1]*phiLy[2]-36.0*rdxCp2[0]*phiLx[2]+((-315.0*rdxCp2[1])-252.0*rdxCp2[0])*phiC[2]+(77.94228634059945*phiC[0]-77.94228634059945*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiLy[3]-80.0*rdxCp2[0]*phiLx[3]+((-105.0*rdxCp2[1])-400.0*rdxCp2[0])*phiC[3]-69.28203230275508*rdxCp2[0]*phiLx[2]-69.28203230275508*rdxCp2[0]*phiC[2]+(25.98076211353316*phiC[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+96.0*rdxCp2[1]*bcVals[3]-10.39230484541326*rdxCp2[1]*phiLy[2]-93.53074360871933*rdxCp2[1]*phiC[2]+((-6.0*phiLy[0])-42.0*phiC[0])*rdxCp2[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+8.660254037844386*rdxCp2[0]*phiC[1]+16.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac-62.35382907247956*rdxCp2[1]*phiLy[3]-561.1844616523159*rdxCp2[1]*phiC[3]+((-36.0*phiLy[1])-252.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiLx[1]-315.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiC[0]-77.94228634059945*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac+8.660254037844386*rdxCp2[0]*phiLx[3]+8.660254037844386*rdxCp2[0]*phiC[3]+138.5640646055102*rdxCp2[1]*bcVals[3]-40.0*rdxCp2[1]*phiLy[2]+9.0*rdxCp2[0]*phiLx[2]+((-200.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+((-34.64101615137754*phiLy[0])-34.64101615137754*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-80.0*rdxCp2[1]*phiLy[3]-25.0*rdxCp2[0]*phiLx[3]+((-400.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+((-69.28203230275508*phiLy[1])-69.28203230275508*phiC[1])*rdxCp2[1]); 

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

  double rdxCp2[2]; 
  double rdxCp2Sq[2]; 
  double rdxCp2R3[2]; 
  double rdxCp2R4[2]; 
  double rdxCp2R6[2]; 
  double rdxCp2R8[2]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxCp2Sq[1]  = rdxCp2[1]*rdxCp2[1]; 
  rdxCp2R3[1]  = rdxCp2[1]*rdxCp2Sq[1]; 
  rdxCp2R4[1]  = rdxCp2Sq[1]*rdxCp2Sq[1]; 
  rdxCp2R6[1]  = rdxCp2Sq[1]*rdxCp2R4[1]; 
  rdxCp2R8[1]  = rdxCp2R4[1]*rdxCp2R4[1]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 
  double *phiUy = phi[3]; 
  double *phiLy = phi[4]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac+16.0*rdxCp2[1]*bcVals[3]+8.660254037844386*rdxCp2[1]*phiLy[2]+8.660254037844386*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]+8.660254037844386*rdxCp2[0]*phiLx[1]+8.660254037844386*rdxCp2[0]*phiC[1]+16.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.01041666666666667*(96.0*rho[1]*volFac+51.96152422706631*rdxCp2[1]*phiLy[3]+51.96152422706631*rdxCp2[1]*phiC[3]+(54.0*phiLy[1]-54.0*phiC[1])*rdxCp2[1]-75.0*rdxCp2[0]*phiLx[1]-315.0*rdxCp2[0]*phiC[1]+138.5640646055102*rdxCp2[0]*bcVals[1]+(77.94228634059945*phiC[0]-77.94228634059945*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.01041666666666667*(96.0*rho[2]*volFac+51.96152422706631*rdxCp2[0]*phiLx[3]+51.96152422706631*rdxCp2[0]*phiC[3]+138.5640646055102*rdxCp2[1]*bcVals[3]-75.0*rdxCp2[1]*phiLy[2]+54.0*rdxCp2[0]*phiLx[2]+((-315.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+(77.94228634059945*phiC[0]-77.94228634059945*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.03125*(32.0*rho[3]*volFac-25.0*rdxCp2[1]*phiLy[3]-25.0*rdxCp2[0]*phiLx[3]+((-105.0*rdxCp2[1])-105.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+(25.98076211353316*phiC[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

}

