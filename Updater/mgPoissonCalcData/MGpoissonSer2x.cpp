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

void MGpoissonGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

void MGpoissonGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((6.928203230275509*rdxCp2[0]*rho[1]+12.0*rho[0]*rdxCp2[1]+100.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+5.196152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]+9.0*phiUx[0]+75.0*phiLy[0]+36.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiUx[1]+(90.0*phiUx[0]+240.0*bcVals[0])*rdxCp2Sq[0])/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[1]+34.64101615137754*rdxCp2[0]*rho[0])*volFac+((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+27.0*rdxCp2[0]*phiLy[1]+(25.98076211353316*phiUy[0]+51.96152422706631*phiUx[0]+25.98076211353316*phiLy[0]-207.8460969082653*bcVals[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiUx[1]+(103.9230484541326*phiUx[0]-207.8460969082653*bcVals[0])*rdxCp2Sq[0])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = ((6.928203230275509*rdxCp2[0]*rho[3]+(80.0*rdxCp2[1]+100.0*rdxCp2[0])*rho[2])*volFac-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[3]+((-69.28203230275508*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-200.0*rdxCp2Sq[1])-250.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(60.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*phiUx[2]+((-200.0*rdxCp2Sq[1])-250.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(173.2050807568877*phiUy[0]-173.2050807568877*phiLy[0])*rdxCp2Sq[1]+(15.0*rdxCp2[0]*phiUy[1]-15.0*rdxCp2[0]*phiLy[1]+(216.5063509461096*phiUy[0]-216.5063509461096*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(800.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[3] = (((160.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[0]*rho[2])*volFac+((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-400.0*rdxCp2[0]*rdxCp2[1])-120.0*rdxCp2Sq[0])*phiUx[3]+((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(346.4101615137754*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiUx[2]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(346.4101615137754*phiUy[1]-346.4101615137754*phiLy[1])*rdxCp2Sq[1]+(77.94228634059945*rdxCp2[0]*phiUy[1]-77.94228634059945*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]-75.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((41.56921938165305*rdxCp2[0]*rho[1]-36.0*rho[0]*rdxCp2[1]-120.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-27.0*phiUy[0])-27.0*phiLy[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+41.56921938165305*rdxCp2[0]*phiUx[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+((-90.0*phiUy[0])-36.0*phiUx[0]-90.0*phiLy[0]+24.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiUx[1]+(160.0*bcVals[0]-60.0*phiUx[0])*rdxCp2Sq[0]))/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[1]-34.64101615137754*rdxCp2[0]*rho[0])*volFac+((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+18.0*rdxCp2[0]*phiLy[1]+((-25.98076211353316*phiUy[0])+51.96152422706631*phiUx[0]-25.98076211353316*phiLy[0]+69.28203230275508*bcVals[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*bcVals[0]*rdxCp2Sq[0])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((20.78460969082652*rdxCp2[0]*rho[3]+((-120.0*rdxCp2[1])-60.0*rdxCp2[0])*rho[2])*volFac-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUy[3]+(138.5640646055102*rdxCp2[0]*rdxCp2[1]+34.64101615137754*rdxCp2Sq[0])*phiUx[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(300.0*rdxCp2Sq[1]+150.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-120.0*rdxCp2[0]*rdxCp2[1])-30.0*rdxCp2Sq[0])*phiUx[2]+(300.0*rdxCp2Sq[1]+150.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(259.8076211353315*phiLy[0]-259.8076211353315*phiUy[0])*rdxCp2Sq[1]+(45.0*rdxCp2[0]*phiUy[1]-45.0*rdxCp2[0]*phiLy[1]+(129.9038105676658*phiLy[0]-129.9038105676658*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(1200.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0]); 
  phiC[3] = (((240.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[0]*rho[2])*volFac+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[2]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiLy[0]-75.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((6.928203230275509*rdxCp2[0]*rho[1]-12.0*rho[0]*rdxCp2[1]-100.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-9.0*phiUy[0])-9.0*phiLy[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]+5.196152422706631*rdxCp2[0]*phiLy[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-36.0*rdxCp2[0]*bcVals[1]+((-75.0*phiUy[0])-75.0*phiLy[0]-9.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiLx[1]-240.0*rdxCp2Sq[0]*bcVals[1]-90.0*phiLx[0]*rdxCp2Sq[0]))/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[1]-34.64101615137754*rdxCp2[0]*rho[0])*volFac+((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]+27.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]+207.8460969082653*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiUy[0])-25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiLx[1]+207.8460969082653*rdxCp2Sq[0]*bcVals[1]-103.9230484541326*phiLx[0]*rdxCp2Sq[0])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((6.928203230275509*rdxCp2[0]*rho[3]+((-80.0*rdxCp2[1])-100.0*rdxCp2[0])*rho[2])*volFac-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-69.28203230275508*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiLx[3]+(200.0*rdxCp2Sq[1]+250.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(200.0*rdxCp2Sq[1]+250.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-60.0*rdxCp2[0]*rdxCp2[1])-90.0*rdxCp2Sq[0])*phiLx[2]+(173.2050807568877*phiLy[0]-173.2050807568877*phiUy[0])*rdxCp2Sq[1]+(15.0*rdxCp2[0]*phiUy[1]-15.0*rdxCp2[0]*phiLy[1]+(216.5063509461096*phiLy[0]-216.5063509461096*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(800.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[3] = (((160.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[0]*rho[2])*volFac+((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-400.0*rdxCp2[0]*rdxCp2[1])-120.0*rdxCp2Sq[0])*phiLx[3]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-346.4101615137754*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*phiUy[1]-346.4101615137754*phiLy[1])*rdxCp2Sq[1]+(77.94228634059945*rdxCp2[0]*phiUy[1]-77.94228634059945*rdxCp2[0]*phiLy[1]+(75.0*phiLy[0]-75.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((41.56921938165305*rdxCp2[0]*rho[1]+36.0*rho[0]*rdxCp2[1]+120.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*phiUy[0]+27.0*phiLy[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+41.56921938165305*rdxCp2[0]*phiLx[1]+24.0*rdxCp2[0]*bcVals[1]+(90.0*phiUy[0]+90.0*phiLy[0]+36.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiLx[1]+160.0*rdxCp2Sq[0]*bcVals[1]+60.0*phiLx[0]*rdxCp2Sq[0])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[1]+34.64101615137754*rdxCp2[0]*rho[0])*volFac+((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]+18.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiUy[0]+25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*bcVals[1])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = ((20.78460969082652*rdxCp2[0]*rho[3]+(120.0*rdxCp2[1]+60.0*rdxCp2[0])*rho[2])*volFac-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUy[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(138.5640646055102*rdxCp2[0]*rdxCp2[1]+34.64101615137754*rdxCp2Sq[0])*phiLx[3]+((-300.0*rdxCp2Sq[1])-150.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-300.0*rdxCp2Sq[1])-150.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(120.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0])*phiLx[2]+(259.8076211353315*phiUy[0]-259.8076211353315*phiLy[0])*rdxCp2Sq[1]+(45.0*rdxCp2[0]*phiUy[1]-45.0*rdxCp2[0]*phiLy[1]+(129.9038105676658*phiUy[0]-129.9038105676658*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1200.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0]); 
  phiC[3] = (((240.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[0]*rho[2])*volFac+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]-346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]-75.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((6.928203230275509*rdxCp2[1]*rho[2]+100.0*rho[0]*rdxCp2[1]+12.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(240.0*rdxCp2Sq[1]+36.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+90.0*phiUy[0]*rdxCp2Sq[1]+((-86.60254037844386*rdxCp2[0]*phiUx[1])+86.60254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUy[0]+75.0*phiUx[0]+75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-10.39230484541326*rdxCp2Sq[0]*phiUx[1]+10.39230484541326*rdxCp2Sq[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0])*rdxCp2Sq[0])/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = ((6.928203230275509*rdxCp2[1]*rho[3]+(100.0*rdxCp2[1]+80.0*rdxCp2[0])*rho[1])*volFac+((-103.9230484541326*rdxCp2Sq[1])-69.28203230275508*rdxCp2[0]*rdxCp2[1])*phiUy[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUx[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLx[3]+15.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-15.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+90.0*phiUy[1]*rdxCp2Sq[1]+(60.0*rdxCp2[0]*phiUy[1]-250.0*rdxCp2[0]*phiUx[1]-250.0*rdxCp2[0]*phiLx[1]+(216.5063509461096*phiUx[0]-216.5063509461096*phiLx[0])*rdxCp2[0])*rdxCp2[1]-200.0*rdxCp2Sq[0]*phiUx[1]-200.0*rdxCp2Sq[0]*phiLx[1]+(173.2050807568877*phiUx[0]-173.2050807568877*phiLx[0])*rdxCp2Sq[0])/(210.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+800.0*rdxCp2Sq[0]); 
  phiC[2] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[2]+34.64101615137754*rho[0]*rdxCp2[1])*volFac+((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]+((-207.8460969082653*rdxCp2Sq[1])-207.8460969082653*rdxCp2[0]*rdxCp2[1])*bcVals[2]+103.9230484541326*phiUy[0]*rdxCp2Sq[1]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]+25.98076211353316*phiUx[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((36.0*rdxCp2[1]+160.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[1]*rho[1])*volFac+((-120.0*rdxCp2Sq[1])-400.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiUx[3]+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiLx[3]+(77.94228634059945*rdxCp2[0]*rdxCp2[1]+346.4101615137754*rdxCp2Sq[0])*phiUx[2]+((-77.94228634059945*rdxCp2[0]*rdxCp2[1])-346.4101615137754*rdxCp2Sq[0])*phiLx[2]+103.9230484541326*phiUy[1]*rdxCp2Sq[1]+(346.4101615137754*rdxCp2[0]*phiUy[1]-86.60254037844386*rdxCp2[0]*phiUx[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiUx[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((41.56921938165305*rdxCp2[1]*rho[2]-120.0*rho[0]*rdxCp2[1]-36.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiUy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]-60.0*phiUy[0]*rdxCp2Sq[1]+(103.9230484541326*rdxCp2[0]*phiUx[1]-103.9230484541326*rdxCp2[0]*phiLx[1]+((-36.0*phiUy[0])-90.0*phiUx[0]-90.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0]*phiUx[1]-31.17691453623978*rdxCp2Sq[0]*phiLx[1]+((-27.0*phiUx[0])-27.0*phiLx[0])*rdxCp2Sq[0]))/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((20.78460969082652*rdxCp2[1]*rho[3]+((-60.0*rdxCp2[1])-120.0*rdxCp2[0])*rho[1])*volFac+(34.64101615137754*rdxCp2Sq[1]+138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiUy[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[3]+45.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-45.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-30.0*phiUy[1]*rdxCp2Sq[1]+((-120.0*rdxCp2[0]*phiUy[1])+150.0*rdxCp2[0]*phiUx[1]+150.0*rdxCp2[0]*phiLx[1]+(129.9038105676658*phiLx[0]-129.9038105676658*phiUx[0])*rdxCp2[0])*rdxCp2[1]+300.0*rdxCp2Sq[0]*phiUx[1]+300.0*rdxCp2Sq[0]*phiLx[1]+(259.8076211353315*phiLx[0]-259.8076211353315*phiUx[0])*rdxCp2Sq[0]))/(30.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+1200.0*rdxCp2Sq[0]); 
  phiC[2] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[2]-34.64101615137754*rho[0]*rdxCp2[1])*volFac+((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]-60.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]-25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((24.0*rdxCp2[1]+240.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[1]*rho[1])*volFac-400.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*rdxCp2[0]*phiUy[1]+86.60254037844386*rdxCp2[0]*phiUx[1]+86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiLx[0]-75.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((6.928203230275509*rdxCp2[1]*rho[2]-100.0*rho[0]*rdxCp2[1]-12.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-240.0*rdxCp2Sq[1])-36.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]-90.0*phiLy[0]*rdxCp2Sq[1]+(86.60254037844386*rdxCp2[0]*phiUx[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+((-75.0*phiUx[0])-9.0*phiLy[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+10.39230484541326*rdxCp2Sq[0]*phiUx[1]-10.39230484541326*rdxCp2Sq[0]*phiLx[1]+((-9.0*phiUx[0])-9.0*phiLx[0])*rdxCp2Sq[0]))/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((6.928203230275509*rdxCp2[1]*rho[3]+((-100.0*rdxCp2[1])-80.0*rdxCp2[0])*rho[1])*volFac-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-103.9230484541326*rdxCp2Sq[1])-69.28203230275508*rdxCp2[0]*rdxCp2[1])*phiLy[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLx[3]+15.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-15.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-90.0*phiLy[1]*rdxCp2Sq[1]+(250.0*rdxCp2[0]*phiUx[1]-60.0*rdxCp2[0]*phiLy[1]+250.0*rdxCp2[0]*phiLx[1]+(216.5063509461096*phiLx[0]-216.5063509461096*phiUx[0])*rdxCp2[0])*rdxCp2[1]+200.0*rdxCp2Sq[0]*phiUx[1]+200.0*rdxCp2Sq[0]*phiLx[1]+(173.2050807568877*phiLx[0]-173.2050807568877*phiUx[0])*rdxCp2Sq[0]))/(210.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+800.0*rdxCp2Sq[0]); 
  phiC[2] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[2]-34.64101615137754*rho[0]*rdxCp2[1])*volFac+((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+(207.8460969082653*rdxCp2Sq[1]+207.8460969082653*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]-103.9230484541326*phiLy[0]*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+((-25.98076211353316*phiUx[0])-51.96152422706631*phiLy[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((36.0*rdxCp2[1]+160.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[1]*rho[1])*volFac+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiUx[3]+((-120.0*rdxCp2Sq[1])-400.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiLx[3]+(77.94228634059945*rdxCp2[0]*rdxCp2[1]+346.4101615137754*rdxCp2Sq[0])*phiUx[2]+((-77.94228634059945*rdxCp2[0]*rdxCp2[1])-346.4101615137754*rdxCp2Sq[0])*phiLx[2]-103.9230484541326*phiLy[1]*rdxCp2Sq[1]+(86.60254037844386*rdxCp2[0]*phiUx[1]-346.4101615137754*rdxCp2[0]*phiLy[1]+86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiLx[0]-75.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((41.56921938165305*rdxCp2[1]*rho[2]+120.0*rho[0]*rdxCp2[1]+36.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiLy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+60.0*phiLy[0]*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2[0]*phiUx[1])+103.9230484541326*rdxCp2[0]*phiLx[1]+(90.0*phiUx[0]+36.0*phiLy[0]+90.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-31.17691453623978*rdxCp2Sq[0]*phiUx[1]+31.17691453623978*rdxCp2Sq[0]*phiLx[1]+(27.0*phiUx[0]+27.0*phiLx[0])*rdxCp2Sq[0])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = ((20.78460969082652*rdxCp2[1]*rho[3]+(60.0*rdxCp2[1]+120.0*rdxCp2[0])*rho[1])*volFac-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[3]+(34.64101615137754*rdxCp2Sq[1]+138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiLy[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[3]+45.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-45.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+30.0*phiLy[1]*rdxCp2Sq[1]+((-150.0*rdxCp2[0]*phiUx[1])+120.0*rdxCp2[0]*phiLy[1]-150.0*rdxCp2[0]*phiLx[1]+(129.9038105676658*phiUx[0]-129.9038105676658*phiLx[0])*rdxCp2[0])*rdxCp2[1]-300.0*rdxCp2Sq[0]*phiUx[1]-300.0*rdxCp2Sq[0]*phiLx[1]+(259.8076211353315*phiUx[0]-259.8076211353315*phiLx[0])*rdxCp2Sq[0])/(30.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+1200.0*rdxCp2Sq[0]); 
  phiC[2] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[2]+34.64101615137754*rho[0]*rdxCp2[1])*volFac+((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]-60.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(25.98076211353316*phiUx[0]-51.96152422706631*phiLy[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((24.0*rdxCp2[1]+240.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[1]*rho[1])*volFac+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+((-86.60254037844386*rdxCp2[0]*phiUx[1])-346.4101615137754*rdxCp2[0]*phiLy[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiUx[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(1454.922678357857*rdxCp2Sq[1]+8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+(8764.17708629852*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-21823.84017536785*rdxCp2R3[1])-134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]-18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(1091.192008768392*rdxCp2[0]*rdxCp2Sq[1]+8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(50400.0*rdxCp2R3[1]+314940.0*rdxCp2[0]*rdxCp2Sq[1]+63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+18900.0*phiUy[0]*rdxCp2R3[1]+(8105.997779422343*rdxCp2[0]*phiUy[1]-18186.53347947321*rdxCp2[0]*phiUx[1]+(116685.0*phiUy[0]+15750.0*phiUx[0]+63000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1091.192008768392*rdxCp2Sq[0]*phiUy[1]-134736.2323207829*rdxCp2Sq[0]*phiUx[1]+(15750.0*phiUy[0]+116685.0*phiUx[0]+314940.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-21823.84017536785*rdxCp2R3[0]*phiUx[1]+(18900.0*phiUx[0]+50400.0*bcVals[0])*rdxCp2R3[0])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = (((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]+3637.306695894642*rdxCp2Sq[0]*rho[0])*volFac+((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-23400.0*rdxCp2[0]*rdxCp2Sq[1])-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(3150.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(50402.67850025433*rdxCp2[0]*rdxCp2Sq[1]+10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+18900.0*phiUy[1]*rdxCp2R3[1]+(20745.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiUx[1]+(20264.99444855586*phiUy[0]+45466.33369868303*phiUx[0]-181865.3347947321*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2835.0*rdxCp2Sq[0]*phiUy[1]-85350.0*rdxCp2Sq[0]*phiUx[1]+(2727.980021920981*phiUy[0]+73915.26821300182*phiUx[0]-164198.4165575295*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-12600.0*rdxCp2R3[0]*phiUx[1]+(10911.92008768392*phiUx[0]-21823.84017536785*bcVals[0])*rdxCp2R3[0])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = (((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+3637.306695894642*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiUx[3]+((-12600.0*rdxCp2R3[1])-85350.0*rdxCp2[0]*rdxCp2Sq[1]-52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(2835.0*rdxCp2[0]*rdxCp2Sq[1]+20745.0*rdxCp2Sq[0]*rdxCp2[1]+18900.0*rdxCp2R3[0])*phiUx[2]+((-21823.84017536785*rdxCp2R3[1])-164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]-181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+10911.92008768392*phiUy[0]*rdxCp2R3[1]+(2160.0*rdxCp2[0]*phiUy[1]-3150.0*rdxCp2[0]*phiUx[1]+(73915.26821300182*phiUy[0]+2727.980021920981*phiUx[0]+10911.92008768392*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(3150.0*rdxCp2Sq[0]*phiUy[1]-23400.0*rdxCp2Sq[0]*phiUx[1]+(45466.33369868303*phiUy[0]+20264.99444855586*phiUx[0]+50402.67850025433*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(320.4293994002423*rdxCp2[0]*rdxCp2[1]+1212.435565298214*rdxCp2Sq[0])*rho[2]+(1212.435565298214*rdxCp2Sq[1]+320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]+1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiUx[3]+((-2078.460969082652*rdxCp2[0]*rdxCp2Sq[1])-3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(2727.980021920981*rdxCp2[0]*rdxCp2Sq[1]+16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]+3637.306695894642*rdxCp2R3[0])*phiUx[2]+(1650.0*rdxCp2[0]*rdxCp2Sq[1]-10500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+3637.306695894642*phiUy[1]*rdxCp2R3[1]+(16116.7327644284*rdxCp2[0]*phiUy[1]-3031.088913245535*rdxCp2[0]*phiUx[1]+(1800.0*phiUy[0]+2625.0*phiUx[0]-10500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2727.980021920981*rdxCp2Sq[0]*phiUy[1]-2078.460969082652*rdxCp2Sq[0]*phiUx[1]+(2625.0*phiUy[0]+1800.0*phiUx[0]+1650.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(415.6921938165305*rdxCp2R3[1]+12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]+26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-2286.307065990918*rdxCp2[0]*rdxCp2Sq[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]-2182.384017536785*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*volFac+(1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(692.8203230275509*rdxCp2R4[1]+21408.14798155132*rdxCp2[0]*rdxCp2R3[1]+61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2R3[1]+11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]+24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(1600.0*rdxCp2R4[1]+48360.0*rdxCp2[0]*rdxCp2R3[1]+111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-600.0*phiUy[0]*rdxCp2R4[1]+((-1039.230484541326*rdxCp2[0]*phiUy[1])+1039.230484541326*rdxCp2[0]*phiUx[1]+((-18540.0*phiUy[0])-900.0*phiUx[0]-3600.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-3065.729929396912*rdxCp2Sq[0]*phiUy[1])+37360.33591926068*rdxCp2Sq[0]*phiUx[1]+((-53070.0*phiUy[0])-32355.0*phiUx[0]-89820.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2182.384017536785*rdxCp2R3[0]*phiUy[1])+94154.28189944413*rdxCp2R3[0]*phiUx[1]+((-31500.0*phiUy[0])-81540.0*phiUx[0]-219960.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+32735.76026305177*rdxCp2R4[0]*phiUx[1]+((-28350.0*phiUx[0])-75600.0*bcVals[0])*rdxCp2R4[0]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]-5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0]*rho[0])*volFac+(692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(3000.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(900.0*rdxCp2[0]*rdxCp2R3[1]+6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6480.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(7967.433714816835*rdxCp2[0]*rdxCp2R3[1]+20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]+3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-600.0*phiUy[1]*rdxCp2R4[1]+((-6240.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiUx[1]+((-2598.076211353316*phiUy[0])-2598.076211353316*phiUx[0]+10392.30484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-16785.0*rdxCp2Sq[0]*phiUy[1])+28650.0*rdxCp2Sq[0]*phiUx[1]+((-7664.324823492281*phiUy[0])-24811.62781842416*phiUx[0]+64951.90528383289*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-5670.0*rdxCp2R3[0]*phiUy[1])+59400.0*rdxCp2R3[0]*phiUx[1]+((-5455.960043841962*phiUy[0])-51441.90898479563*phiUx[0]+113795.7380572752*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+18900.0*rdxCp2R4[0]*phiUx[1]+(32735.76026305177*bcVals[0]-16367.88013152588*phiUx[0])*rdxCp2R4[0]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = (((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+(240.0*rdxCp2R3[1]+7740.0*rdxCp2[0]*rdxCp2Sq[1]+30300.0*rdxCp2Sq[0]*rdxCp2[1]+31500.0*rdxCp2R3[0])*rho[2]+((-720.0*rdxCp2[0]*rdxCp2Sq[1])-1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiUx[3]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-52500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(180.0*rdxCp2[0]*rdxCp2R3[1]+6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25380.0*rdxCp2R3[0]*rdxCp2[1]+28350.0*rdxCp2R4[0])*phiUx[2]+(692.8203230275509*rdxCp2R4[1]+21823.84017536785*rdxCp2[0]*rdxCp2R3[1]+72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]+60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-300.0*rdxCp2[0]*phiUy[1])+300.0*rdxCp2[0]*phiUx[1]+(779.4228634059946*phiUy[0]-259.8076211353315*phiUx[0]-1039.230484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(10800.0*rdxCp2Sq[0]*phiUx[1]+(21823.84017536785*phiUy[0]-9353.074360871933*phiUx[0]-24941.53162899183*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3150.0*rdxCp2R3[0]*phiUy[1]+23400.0*rdxCp2R3[0]*phiUx[1]+(45466.33369868303*phiUy[0]-20264.99444855586*phiUx[0]-50402.67850025433*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiUx[3]+(866.0254037844386*rdxCp2[0]*rdxCp2R3[1]-9093.266739736604*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(519.6152422706631*rdxCp2[0]*rdxCp2R3[1]+9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]+36476.99000740053*rdxCp2R3[0]*rdxCp2[1]+16367.88013152588*rdxCp2R4[0])*phiUx[2]+(2600.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(4330.127018922193*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiUx[1]+((-750.0*phiUy[0])-750.0*phiUx[0]+3000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]+6235.382907247957*rdxCp2Sq[0]*phiUx[1]+(10800.0*bcVals[0]-5400.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(8183.940065762942*rdxCp2R3[0]*phiUy[1]+6235.382907247957*rdxCp2R3[0]*phiUx[1]+(7875.0*phiUy[0]-5400.0*phiUx[0]-4950.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2182.384017536785*rdxCp2R3[1])-6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]-2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(26292.53125889555*rdxCp2[0]*rdxCp2Sq[1]+12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*volFac+((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(32735.76026305177*rdxCp2R4[1]+94154.28189944413*rdxCp2[0]*rdxCp2R3[1]+37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-2182.384017536785*rdxCp2[0]*rdxCp2R3[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-75600.0*rdxCp2R4[1])-219960.0*rdxCp2[0]*rdxCp2R3[1]-89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-28350.0*phiUy[0]*rdxCp2R4[1]+(24317.99333826703*rdxCp2[0]*phiUy[1]+36373.06695894642*rdxCp2[0]*phiUx[1]+((-81540.0*phiUy[0])-31500.0*phiUx[0]+21000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(11223.68923304632*rdxCp2Sq[0]*phiUy[1]+61279.95757178687*rdxCp2Sq[0]*phiUx[1]+((-32355.0*phiUy[0])-53070.0*phiUx[0]+111280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2R3[0]*phiUy[1]+21408.14798155132*rdxCp2R3[0]*phiUx[1]+((-900.0*phiUy[0])-18540.0*phiUx[0]+48360.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+692.8203230275509*rdxCp2R4[0]*phiUx[1]+(1600.0*bcVals[0]-600.0*phiUx[0])*rdxCp2R4[0]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = (((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1770.0*rdxCp2[0]*rdxCp2Sq[1])-720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(31500.0*rdxCp2R3[1]+30300.0*rdxCp2[0]*rdxCp2Sq[1]+7740.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(23400.0*rdxCp2[0]*rdxCp2R3[1]+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3150.0*rdxCp2[0]*rdxCp2R3[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-50402.67850025433*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+28350.0*phiUy[1]*rdxCp2R4[1]+(25380.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiUx[1]+((-20264.99444855586*phiUy[0])+45466.33369868303*phiUx[0]+60621.7782649107*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(6435.0*rdxCp2Sq[0]*phiUy[1]-25200.0*rdxCp2Sq[0]*phiUx[1]+((-9353.074360871933*phiUy[0])+21823.84017536785*phiUx[0]+72919.33899864974*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(180.0*rdxCp2R3[0]*phiUy[1]-900.0*rdxCp2R3[0]*phiUx[1]+((-259.8076211353315*phiUy[0])+779.4228634059946*phiUx[0]+21823.84017536785*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+692.8203230275509*bcVals[0]*rdxCp2R4[0])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-5455.960043841962*rho[0]*rdxCp2R3[1]-15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiUx[3]+(18900.0*rdxCp2R4[1]+59400.0*rdxCp2[0]*rdxCp2R3[1]+28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-5670.0*rdxCp2[0]*rdxCp2R3[1])-16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiUx[2]+(32735.76026305177*rdxCp2R4[1]+113795.7380572752*rdxCp2[0]*rdxCp2R3[1]+64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-16367.88013152588*phiUy[0]*rdxCp2R4[1]+(6480.0*rdxCp2[0]*phiUy[1]+6300.0*rdxCp2[0]*phiUx[1]+((-51441.90898479563*phiUy[0])-5455.960043841962*phiUx[0]+3637.306695894642*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(6480.0*rdxCp2Sq[0]*phiUy[1]+8850.0*rdxCp2Sq[0]*phiUx[1]+((-24811.62781842416*phiUy[0])-7664.324823492281*phiUx[0]+20438.19952931275*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(900.0*rdxCp2R3[0]*phiUy[1]+3000.0*rdxCp2R3[0]*phiUx[1]+((-2598.076211353316*phiUy[0])-2598.076211353316*phiUx[0]+7967.433714816835*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(6235.382907247957*rdxCp2[0]*rdxCp2R3[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(8183.940065762942*rdxCp2[0]*rdxCp2R3[1]+21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]+4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-4950.0*rdxCp2[0]*rdxCp2R3[1])+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+16367.88013152588*phiUy[1]*rdxCp2R4[1]+(36476.99000740053*rdxCp2[0]*phiUy[1]-9093.266739736604*rdxCp2[0]*phiUx[1]+((-5400.0*phiUy[0])+7875.0*phiUx[0]+10500.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(9846.708841029065*rdxCp2Sq[0]*phiUy[1]+(8850.0*bcVals[0]-5400.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2R3[0]*phiUy[1]+866.0254037844386*rdxCp2R3[0]*phiUx[1]+((-750.0*phiUy[0])-750.0*phiUx[0]+2600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-51.96152422706631*rdxCp2Sq[1])-571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+((-571.5767664977294*rdxCp2[0]*rdxCp2[1])-51.96152422706631*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*volFac+(300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-86.60254037844386*rdxCp2R3[1])-987.26896031426*rdxCp2[0]*rdxCp2Sq[1]-173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2Sq[1])-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-200.0*rdxCp2R3[1])-2220.0*rdxCp2[0]*rdxCp2Sq[1]-100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+75.0*phiUy[0]*rdxCp2R3[1]+((-259.8076211353315*rdxCp2[0]*phiUy[1])-173.2050807568877*rdxCp2[0]*phiUx[1]+(855.0*phiUy[0]+150.0*phiUx[0]-100.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-51.96152422706631*rdxCp2Sq[0]*phiUy[1])-987.26896031426*rdxCp2Sq[0]*phiUx[1]+(150.0*phiUy[0]+855.0*phiUx[0]-2220.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-86.60254037844386*rdxCp2R3[0]*phiUx[1]+(75.0*phiUx[0]-200.0*bcVals[0])*rdxCp2R3[0])/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0]*rho[0])*volFac+(259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(225.0*rdxCp2[0]*rdxCp2Sq[1]-225.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1991.858428704209*rdxCp2[0]*rdxCp2Sq[1])-86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-225.0*phiUy[1]*rdxCp2R3[1]+((-765.0*rdxCp2[0]*phiUy[1])+750.0*rdxCp2[0]*phiUx[1]+(649.5190528383289*phiUy[0]-649.5190528383289*phiUx[0]-866.0254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90.0*rdxCp2Sq[0]*phiUy[1])+150.0*rdxCp2Sq[0]*phiUx[1]+(129.9038105676658*phiUy[0]-129.9038105676658*phiUx[0]-3031.088913245535*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-259.8076211353315*bcVals[0]*rdxCp2R3[0]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+129.9038105676658*rho[0]*rdxCp2Sq[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiUx[3]+(150.0*rdxCp2[0]*rdxCp2Sq[1]+750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-90.0*rdxCp2[0]*rdxCp2Sq[1])-765.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiUx[2]+((-259.8076211353315*rdxCp2R3[1])-3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-225.0*rdxCp2[0]*phiUy[1])-150.0*rdxCp2[0]*phiUx[1]+((-129.9038105676658*phiUy[0])+129.9038105676658*phiUx[0]-86.60254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(225.0*rdxCp2Sq[0]*phiUy[1]-750.0*rdxCp2Sq[0]*phiUx[1]+((-649.5190528383289*phiUy[0])+649.5190528383289*phiUx[0]-1991.858428704209*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+((-363.7306695894642*rdxCp2[0]*rdxCp2[1])-259.8076211353315*rdxCp2Sq[0])*rho[2]+((-259.8076211353315*rdxCp2Sq[1])-363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]+1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1300.0*rdxCp2[0]*rdxCp2Sq[1])-500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1299.038105676658*rdxCp2[0]*phiUy[1]+433.0127018922193*rdxCp2[0]*phiUx[1]+(375.0*phiUy[0]-375.0*phiUx[0]-500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(259.8076211353315*rdxCp2Sq[0]*phiUy[1]-433.0127018922193*rdxCp2Sq[0]*phiUx[1]+((-375.0*phiUy[0])+375.0*phiUx[0]-1300.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(1454.922678357857*rdxCp2Sq[1]+8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+((-8764.17708629852*rdxCp2[0]*rdxCp2[1])-1454.922678357857*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-50400.0*rdxCp2R3[1])-314940.0*rdxCp2[0]*rdxCp2Sq[1]-63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1091.192008768392*rdxCp2[0]*rdxCp2Sq[1]+8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-21823.84017536785*rdxCp2R3[1])-134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]-18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-18900.0*phiLy[0]*rdxCp2R3[1]+(18186.53347947321*rdxCp2[0]*phiUx[1]-8105.997779422343*rdxCp2[0]*phiLy[1]+((-15750.0*phiUx[0])-116685.0*phiLy[0]-63000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(134736.2323207829*rdxCp2Sq[0]*phiUx[1]-1091.192008768392*rdxCp2Sq[0]*phiLy[1]+((-116685.0*phiUx[0])-15750.0*phiLy[0]-314940.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+21823.84017536785*rdxCp2R3[0]*phiUx[1]+((-18900.0*phiUx[0])-50400.0*bcVals[0])*rdxCp2R3[0]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]-3637.306695894642*rdxCp2Sq[0]*rho[0])*volFac+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-50402.67850025433*rdxCp2[0]*rdxCp2Sq[1])-10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(3150.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-23400.0*rdxCp2[0]*rdxCp2Sq[1])-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-18900.0*phiLy[1]*rdxCp2R3[1]+(52500.0*rdxCp2[0]*phiUx[1]-20745.0*rdxCp2[0]*phiLy[1]+((-45466.33369868303*phiUx[0])-20264.99444855586*phiLy[0]+181865.3347947321*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(85350.0*rdxCp2Sq[0]*phiUx[1]-2835.0*rdxCp2Sq[0]*phiLy[1]+((-73915.26821300182*phiUx[0])-2727.980021920981*phiLy[0]+164198.4165575295*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+12600.0*rdxCp2R3[0]*phiUx[1]+(21823.84017536785*bcVals[0]-10911.92008768392*phiUx[0])*rdxCp2R3[0]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = (((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-3637.306695894642*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiUx[3]+((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(21823.84017536785*rdxCp2R3[1]+164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]+181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2835.0*rdxCp2[0]*rdxCp2Sq[1]+20745.0*rdxCp2Sq[0]*rdxCp2[1]+18900.0*rdxCp2R3[0])*phiUx[2]+((-12600.0*rdxCp2R3[1])-85350.0*rdxCp2[0]*rdxCp2Sq[1]-52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-10911.92008768392*phiLy[0]*rdxCp2R3[1]+(3150.0*rdxCp2[0]*phiUx[1]-2160.0*rdxCp2[0]*phiLy[1]+((-2727.980021920981*phiUx[0])-73915.26821300182*phiLy[0]-10911.92008768392*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(23400.0*rdxCp2Sq[0]*phiUx[1]-3150.0*rdxCp2Sq[0]*phiLy[1]+((-20264.99444855586*phiUx[0])-45466.33369868303*phiLy[0]-50402.67850025433*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(320.4293994002423*rdxCp2[0]*rdxCp2[1]+1212.435565298214*rdxCp2Sq[0])*rho[2]+((-1212.435565298214*rdxCp2Sq[1])-320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]-1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiUx[3]+((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(10500.0*rdxCp2Sq[0]*rdxCp2[1]-1650.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(2727.980021920981*rdxCp2[0]*rdxCp2Sq[1]+16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]+3637.306695894642*rdxCp2R3[0])*phiUx[2]+((-2078.460969082652*rdxCp2[0]*rdxCp2Sq[1])-3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-3637.306695894642*phiLy[1]*rdxCp2R3[1]+(3031.088913245535*rdxCp2[0]*phiUx[1]-16116.7327644284*rdxCp2[0]*phiLy[1]+((-2625.0*phiUx[0])-1800.0*phiLy[0]+10500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2078.460969082652*rdxCp2Sq[0]*phiUx[1]-2727.980021920981*rdxCp2Sq[0]*phiLy[1]+((-1800.0*phiUx[0])-2625.0*phiLy[0]-1650.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(415.6921938165305*rdxCp2R3[1]+12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]+26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(2286.307065990918*rdxCp2[0]*rdxCp2Sq[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*volFac+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(1600.0*rdxCp2R4[1]+48360.0*rdxCp2[0]*rdxCp2R3[1]+111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(311.7691453623978*rdxCp2[0]*rdxCp2R3[1]+11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]+24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(692.8203230275509*rdxCp2R4[1]+21408.14798155132*rdxCp2[0]*rdxCp2R3[1]+61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+600.0*phiLy[0]*rdxCp2R4[1]+((-1039.230484541326*rdxCp2[0]*phiUx[1])+1039.230484541326*rdxCp2[0]*phiLy[1]+(900.0*phiUx[0]+18540.0*phiLy[0]+3600.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-37360.33591926068*rdxCp2Sq[0]*phiUx[1])+3065.729929396912*rdxCp2Sq[0]*phiLy[1]+(32355.0*phiUx[0]+53070.0*phiLy[0]+89820.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-94154.28189944413*rdxCp2R3[0]*phiUx[1])+2182.384017536785*rdxCp2R3[0]*phiLy[1]+(81540.0*phiUx[0]+31500.0*phiLy[0]+219960.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-32735.76026305177*rdxCp2R4[0]*phiUx[1]+(28350.0*phiUx[0]+75600.0*bcVals[0])*rdxCp2R4[0])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = (((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]+5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0]*rho[0])*volFac+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(7967.433714816835*rdxCp2[0]*rdxCp2R3[1]+20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]+3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(900.0*rdxCp2[0]*rdxCp2R3[1]+6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6480.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(3000.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+600.0*phiLy[1]*rdxCp2R4[1]+((-3000.0*rdxCp2[0]*phiUx[1])+6240.0*rdxCp2[0]*phiLy[1]+(2598.076211353316*phiUx[0]+2598.076211353316*phiLy[0]-10392.30484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-28650.0*rdxCp2Sq[0]*phiUx[1])+16785.0*rdxCp2Sq[0]*phiLy[1]+(24811.62781842416*phiUx[0]+7664.324823492281*phiLy[0]-64951.90528383289*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-59400.0*rdxCp2R3[0]*phiUx[1])+5670.0*rdxCp2R3[0]*phiLy[1]+(51441.90898479563*phiUx[0]+5455.960043841962*phiLy[0]-113795.7380572752*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-18900.0*rdxCp2R4[0]*phiUx[1]+(16367.88013152588*phiUx[0]-32735.76026305177*bcVals[0])*rdxCp2R4[0])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = (((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+(240.0*rdxCp2R3[1]+7740.0*rdxCp2[0]*rdxCp2Sq[1]+30300.0*rdxCp2Sq[0]*rdxCp2[1]+31500.0*rdxCp2R3[0])*rho[2]+(720.0*rdxCp2[0]*rdxCp2Sq[1]+1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiUx[3]+(346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(692.8203230275509*rdxCp2R4[1]+21823.84017536785*rdxCp2[0]*rdxCp2R3[1]+72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]+60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(180.0*rdxCp2[0]*rdxCp2R3[1]+6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25380.0*rdxCp2R3[0]*rdxCp2[1]+28350.0*rdxCp2R4[0])*phiUx[2]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-52500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-300.0*rdxCp2[0]*phiUx[1])+300.0*rdxCp2[0]*phiLy[1]+(259.8076211353315*phiUx[0]-779.4228634059946*phiLy[0]+1039.230484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((9353.074360871933*phiUx[0]-21823.84017536785*phiLy[0]+24941.53162899183*bcVals[0])*rdxCp2Sq[0]-10800.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-23400.0*rdxCp2R3[0]*phiUx[1])-3150.0*rdxCp2R3[0]*phiLy[1]+(20264.99444855586*phiUx[0]-45466.33369868303*phiLy[0]+50402.67850025433*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiUx[3]+((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(2600.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(519.6152422706631*rdxCp2[0]*rdxCp2R3[1]+9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]+36476.99000740053*rdxCp2R3[0]*rdxCp2[1]+16367.88013152588*rdxCp2R4[0])*phiUx[2]+(866.0254037844386*rdxCp2[0]*rdxCp2R3[1]-9093.266739736604*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-866.0254037844386*rdxCp2[0]*phiUx[1])-4330.127018922193*rdxCp2[0]*phiLy[1]+(750.0*phiUx[0]+750.0*phiLy[0]-3000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-6235.382907247957*rdxCp2Sq[0]*phiUx[1])-21823.84017536785*rdxCp2Sq[0]*phiLy[1]+(5400.0*phiUx[0]-10800.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-6235.382907247957*rdxCp2R3[0]*phiUx[1])-8183.940065762942*rdxCp2R3[0]*phiLy[1]+(5400.0*phiUx[0]-7875.0*phiLy[0]+4950.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2182.384017536785*rdxCp2R3[1])-6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]-2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-26292.53125889555*rdxCp2[0]*rdxCp2Sq[1])-12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]-415.6921938165305*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*volFac+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(75600.0*rdxCp2R4[1]+219960.0*rdxCp2[0]*rdxCp2R3[1]+89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-2182.384017536785*rdxCp2[0]*rdxCp2R3[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(32735.76026305177*rdxCp2R4[1]+94154.28189944413*rdxCp2[0]*rdxCp2R3[1]+37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+28350.0*phiLy[0]*rdxCp2R4[1]+((-36373.06695894642*rdxCp2[0]*phiUx[1])-24317.99333826703*rdxCp2[0]*phiLy[1]+(31500.0*phiUx[0]+81540.0*phiLy[0]-21000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-61279.95757178687*rdxCp2Sq[0]*phiUx[1])-11223.68923304632*rdxCp2Sq[0]*phiLy[1]+(53070.0*phiUx[0]+32355.0*phiLy[0]-111280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-21408.14798155132*rdxCp2R3[0]*phiUx[1])-311.7691453623978*rdxCp2R3[0]*phiLy[1]+(18540.0*phiUx[0]+900.0*phiLy[0]-48360.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-692.8203230275509*rdxCp2R4[0]*phiUx[1]+(600.0*phiUx[0]-1600.0*bcVals[0])*rdxCp2R4[0])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1770.0*rdxCp2[0]*rdxCp2Sq[1])-720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-31500.0*rdxCp2R3[1])-30300.0*rdxCp2[0]*rdxCp2Sq[1]-7740.0*rdxCp2Sq[0]*rdxCp2[1]-240.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(50402.67850025433*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3150.0*rdxCp2[0]*rdxCp2R3[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(23400.0*rdxCp2[0]*rdxCp2R3[1]+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-28350.0*phiLy[1]*rdxCp2R4[1]+(52500.0*rdxCp2[0]*phiUx[1]-25380.0*rdxCp2[0]*phiLy[1]+((-45466.33369868303*phiUx[0])+20264.99444855586*phiLy[0]-60621.7782649107*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(25200.0*rdxCp2Sq[0]*phiUx[1]-6435.0*rdxCp2Sq[0]*phiLy[1]+((-21823.84017536785*phiUx[0])+9353.074360871933*phiLy[0]-72919.33899864974*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(900.0*rdxCp2R3[0]*phiUx[1]-180.0*rdxCp2R3[0]*phiLy[1]+((-779.4228634059946*phiUx[0])+259.8076211353315*phiLy[0]-21823.84017536785*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-692.8203230275509*bcVals[0]*rdxCp2R4[0]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+5455.960043841962*rho[0]*rdxCp2R3[1]+15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiUx[3]+((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-32735.76026305177*rdxCp2R4[1])-113795.7380572752*rdxCp2[0]*rdxCp2R3[1]-64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-5670.0*rdxCp2[0]*rdxCp2R3[1])-16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiUx[2]+(18900.0*rdxCp2R4[1]+59400.0*rdxCp2[0]*rdxCp2R3[1]+28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+16367.88013152588*phiLy[0]*rdxCp2R4[1]+((-6300.0*rdxCp2[0]*phiUx[1])-6480.0*rdxCp2[0]*phiLy[1]+(5455.960043841962*phiUx[0]+51441.90898479563*phiLy[0]-3637.306695894642*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-8850.0*rdxCp2Sq[0]*phiUx[1])-6480.0*rdxCp2Sq[0]*phiLy[1]+(7664.324823492281*phiUx[0]+24811.62781842416*phiLy[0]-20438.19952931275*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3000.0*rdxCp2R3[0]*phiUx[1])-900.0*rdxCp2R3[0]*phiLy[1]+(2598.076211353316*phiUx[0]+2598.076211353316*phiLy[0]-7967.433714816835*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(4950.0*rdxCp2[0]*rdxCp2R3[1]-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(8183.940065762942*rdxCp2[0]*rdxCp2R3[1]+21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]+4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(6235.382907247957*rdxCp2[0]*rdxCp2R3[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-16367.88013152588*phiLy[1]*rdxCp2R4[1]+(9093.266739736604*rdxCp2[0]*phiUx[1]-36476.99000740053*rdxCp2[0]*phiLy[1]+((-7875.0*phiUx[0])+5400.0*phiLy[0]-10500.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((5400.0*phiLy[0]-8850.0*bcVals[0])*rdxCp2Sq[0]-9846.708841029065*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-866.0254037844386*rdxCp2R3[0]*phiUx[1])-519.6152422706631*rdxCp2R3[0]*phiLy[1]+(750.0*phiUx[0]+750.0*phiLy[0]-2600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-51.96152422706631*rdxCp2Sq[1])-571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+(571.5767664977294*rdxCp2[0]*rdxCp2[1]+51.96152422706631*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*volFac+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-200.0*rdxCp2R3[1])-2220.0*rdxCp2[0]*rdxCp2Sq[1]-100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-51.96152422706631*rdxCp2[0]*rdxCp2Sq[1])-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-86.60254037844386*rdxCp2R3[1])-987.26896031426*rdxCp2[0]*rdxCp2Sq[1]-173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-75.0*phiLy[0]*rdxCp2R3[1]+(173.2050807568877*rdxCp2[0]*phiUx[1]+259.8076211353315*rdxCp2[0]*phiLy[1]+((-150.0*phiUx[0])-855.0*phiLy[0]+100.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(987.26896031426*rdxCp2Sq[0]*phiUx[1]+51.96152422706631*rdxCp2Sq[0]*phiLy[1]+((-855.0*phiUx[0])-150.0*phiLy[0]+2220.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+86.60254037844386*rdxCp2R3[0]*phiUx[1]+(200.0*bcVals[0]-75.0*phiUx[0])*rdxCp2R3[0]))/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = (((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]-129.9038105676658*rdxCp2Sq[0]*rho[0])*volFac+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1991.858428704209*rdxCp2[0]*rdxCp2Sq[1])-86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(225.0*rdxCp2[0]*rdxCp2Sq[1]-225.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+225.0*phiLy[1]*rdxCp2R3[1]+((-750.0*rdxCp2[0]*phiUx[1])+765.0*rdxCp2[0]*phiLy[1]+(649.5190528383289*phiUx[0]-649.5190528383289*phiLy[0]+866.0254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-150.0*rdxCp2Sq[0]*phiUx[1])+90.0*rdxCp2Sq[0]*phiLy[1]+(129.9038105676658*phiUx[0]-129.9038105676658*phiLy[0]+3031.088913245535*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+259.8076211353315*bcVals[0]*rdxCp2R3[0])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-129.9038105676658*rho[0]*rdxCp2Sq[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiUx[3]+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-259.8076211353315*rdxCp2R3[1])-3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-90.0*rdxCp2[0]*rdxCp2Sq[1])-765.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiUx[2]+(150.0*rdxCp2[0]*rdxCp2Sq[1]+750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(150.0*rdxCp2[0]*phiUx[1]+225.0*rdxCp2[0]*phiLy[1]+((-129.9038105676658*phiUx[0])+129.9038105676658*phiLy[0]+86.60254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(750.0*rdxCp2Sq[0]*phiUx[1]-225.0*rdxCp2Sq[0]*phiLy[1]+((-649.5190528383289*phiUx[0])+649.5190528383289*phiLy[0]+1991.858428704209*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+((-363.7306695894642*rdxCp2[0]*rdxCp2[1])-259.8076211353315*rdxCp2Sq[0])*rho[2]+(259.8076211353315*rdxCp2Sq[1]+363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1300.0*rdxCp2[0]*rdxCp2Sq[1])-500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]+1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-433.0127018922193*rdxCp2[0]*phiUx[1])-1299.038105676658*rdxCp2[0]*phiLy[1]+(375.0*phiUx[0]-375.0*phiLy[0]+500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(433.0127018922193*rdxCp2Sq[0]*phiUx[1]-259.8076211353315*rdxCp2Sq[0]*phiLy[1]+((-375.0*phiUx[0])+375.0*phiLy[0]+1300.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-1454.922678357857*rdxCp2Sq[1])-8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+(8764.17708629852*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(21823.84017536785*rdxCp2R3[1]+134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]+18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-1091.192008768392*rdxCp2[0]*rdxCp2Sq[1])-8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-50400.0*rdxCp2R3[1])-314940.0*rdxCp2[0]*rdxCp2Sq[1]-63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-18900.0*phiUy[0]*rdxCp2R3[1]+(8105.997779422343*rdxCp2[0]*phiUy[1]-18186.53347947321*rdxCp2[0]*phiLx[1]-63000.0*rdxCp2[0]*bcVals[1]+((-116685.0*phiUy[0])-15750.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(1091.192008768392*rdxCp2Sq[0]*phiUy[1]-134736.2323207829*rdxCp2Sq[0]*phiLx[1]-314940.0*rdxCp2Sq[0]*bcVals[1]+((-15750.0*phiUy[0])-116685.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-21823.84017536785*rdxCp2R3[0]*phiLx[1]-50400.0*rdxCp2R3[0]*bcVals[1]-18900.0*phiLx[0]*rdxCp2R3[0]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = (((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]-3637.306695894642*rdxCp2Sq[0]*rho[0])*volFac+((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(23400.0*rdxCp2[0]*rdxCp2Sq[1]+3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-50402.67850025433*rdxCp2[0]*rdxCp2Sq[1])-10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+18900.0*phiUy[1]*rdxCp2R3[1]+(20745.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiLx[1]+181865.3347947321*rdxCp2[0]*bcVals[1]+((-20264.99444855586*phiUy[0])-45466.33369868303*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(2835.0*rdxCp2Sq[0]*phiUy[1]-85350.0*rdxCp2Sq[0]*phiLx[1]+164198.4165575295*rdxCp2Sq[0]*bcVals[1]+((-2727.980021920981*phiUy[0])-73915.26821300182*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-12600.0*rdxCp2R3[0]*phiLx[1]+21823.84017536785*rdxCp2R3[0]*bcVals[1]-10911.92008768392*phiLx[0]*rdxCp2R3[0])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-3637.306695894642*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiLx[3]+(12600.0*rdxCp2R3[1]+85350.0*rdxCp2[0]*rdxCp2Sq[1]+52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-2835.0*rdxCp2[0]*rdxCp2Sq[1])-20745.0*rdxCp2Sq[0]*rdxCp2[1]-18900.0*rdxCp2R3[0])*phiLx[2]+(21823.84017536785*rdxCp2R3[1]+164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]+181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-10911.92008768392*phiUy[0]*rdxCp2R3[1]+(2160.0*rdxCp2[0]*phiUy[1]-3150.0*rdxCp2[0]*phiLx[1]-10911.92008768392*rdxCp2[0]*bcVals[1]+((-73915.26821300182*phiUy[0])-2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(3150.0*rdxCp2Sq[0]*phiUy[1]-23400.0*rdxCp2Sq[0]*phiLx[1]-50402.67850025433*rdxCp2Sq[0]*bcVals[1]+((-45466.33369868303*phiUy[0])-20264.99444855586*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-320.4293994002423*rdxCp2[0]*rdxCp2[1])-1212.435565298214*rdxCp2Sq[0])*rho[2]+(1212.435565298214*rdxCp2Sq[1]+320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]-1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiLx[3]+(2078.460969082652*rdxCp2[0]*rdxCp2Sq[1]+3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-2727.980021920981*rdxCp2[0]*rdxCp2Sq[1])-16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]-3637.306695894642*rdxCp2R3[0])*phiLx[2]+(10500.0*rdxCp2Sq[0]*rdxCp2[1]-1650.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+3637.306695894642*phiUy[1]*rdxCp2R3[1]+(16116.7327644284*rdxCp2[0]*phiUy[1]-3031.088913245535*rdxCp2[0]*phiLx[1]+10500.0*rdxCp2[0]*bcVals[1]+((-1800.0*phiUy[0])-2625.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(2727.980021920981*rdxCp2Sq[0]*phiUy[1]-2078.460969082652*rdxCp2Sq[0]*phiLx[1]-1650.0*rdxCp2Sq[0]*bcVals[1]+((-2625.0*phiUy[0])-1800.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-415.6921938165305*rdxCp2R3[1])-12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]-26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-2286.307065990918*rdxCp2[0]*rdxCp2Sq[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]-2182.384017536785*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*volFac+(1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-692.8203230275509*rdxCp2R4[1])-21408.14798155132*rdxCp2[0]*rdxCp2R3[1]-61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2R3[1])-11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]-24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-1600.0*rdxCp2R4[1])-48360.0*rdxCp2[0]*rdxCp2R3[1]-111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+600.0*phiUy[0]*rdxCp2R4[1]+((-1039.230484541326*rdxCp2[0]*phiUy[1])+1039.230484541326*rdxCp2[0]*phiLx[1]+3600.0*rdxCp2[0]*bcVals[1]+(18540.0*phiUy[0]+900.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-3065.729929396912*rdxCp2Sq[0]*phiUy[1])+37360.33591926068*rdxCp2Sq[0]*phiLx[1]+89820.0*rdxCp2Sq[0]*bcVals[1]+(53070.0*phiUy[0]+32355.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2182.384017536785*rdxCp2R3[0]*phiUy[1])+94154.28189944413*rdxCp2R3[0]*phiLx[1]+219960.0*rdxCp2R3[0]*bcVals[1]+(31500.0*phiUy[0]+81540.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+32735.76026305177*rdxCp2R4[0]*phiLx[1]+75600.0*rdxCp2R4[0]*bcVals[1]+28350.0*phiLx[0]*rdxCp2R4[0])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]+5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0]*rho[0])*volFac+(692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3000.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6480.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-7967.433714816835*rdxCp2[0]*rdxCp2R3[1])-20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-600.0*phiUy[1]*rdxCp2R4[1]+((-6240.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiLx[1]-10392.30484541326*rdxCp2[0]*bcVals[1]+(2598.076211353316*phiUy[0]+2598.076211353316*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-16785.0*rdxCp2Sq[0]*phiUy[1])+28650.0*rdxCp2Sq[0]*phiLx[1]-64951.90528383289*rdxCp2Sq[0]*bcVals[1]+(7664.324823492281*phiUy[0]+24811.62781842416*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-5670.0*rdxCp2R3[0]*phiUy[1])+59400.0*rdxCp2R3[0]*phiLx[1]-113795.7380572752*rdxCp2R3[0]*bcVals[1]+(5455.960043841962*phiUy[0]+51441.90898479563*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+18900.0*rdxCp2R4[0]*phiLx[1]-32735.76026305177*rdxCp2R4[0]*bcVals[1]+16367.88013152588*phiLx[0]*rdxCp2R4[0]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+((-240.0*rdxCp2R3[1])-7740.0*rdxCp2[0]*rdxCp2Sq[1]-30300.0*rdxCp2Sq[0]*rdxCp2[1]-31500.0*rdxCp2R3[0])*rho[2]+((-720.0*rdxCp2[0]*rdxCp2Sq[1])-1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiLx[3]+(900.0*rdxCp2[0]*rdxCp2R3[1]+25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+52500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-180.0*rdxCp2[0]*rdxCp2R3[1])-6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25380.0*rdxCp2R3[0]*rdxCp2[1]-28350.0*rdxCp2R4[0])*phiLx[2]+((-692.8203230275509*rdxCp2R4[1])-21823.84017536785*rdxCp2[0]*rdxCp2R3[1]-72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]-60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-300.0*rdxCp2[0]*phiUy[1])+300.0*rdxCp2[0]*phiLx[1]+1039.230484541326*rdxCp2[0]*bcVals[1]+(259.8076211353315*phiLx[0]-779.4228634059946*phiUy[0])*rdxCp2[0])*rdxCp2R3[1]+(10800.0*rdxCp2Sq[0]*phiLx[1]+24941.53162899183*rdxCp2Sq[0]*bcVals[1]+(9353.074360871933*phiLx[0]-21823.84017536785*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3150.0*rdxCp2R3[0]*phiUy[1]+23400.0*rdxCp2R3[0]*phiLx[1]+50402.67850025433*rdxCp2R3[0]*bcVals[1]+(20264.99444855586*phiLx[0]-45466.33369868303*phiUy[0])*rdxCp2R3[0])*rdxCp2[1]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiLx[3]+(9093.266739736604*rdxCp2R3[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2R3[1])-9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]-36476.99000740053*rdxCp2R3[0]*rdxCp2[1]-16367.88013152588*rdxCp2R4[0])*phiLx[2]+((-2600.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(4330.127018922193*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiLx[1]-3000.0*rdxCp2[0]*bcVals[1]+(750.0*phiUy[0]+750.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]+6235.382907247957*rdxCp2Sq[0]*phiLx[1]-10800.0*rdxCp2Sq[0]*bcVals[1]+5400.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(8183.940065762942*rdxCp2R3[0]*phiUy[1]+6235.382907247957*rdxCp2R3[0]*phiLx[1]+4950.0*rdxCp2R3[0]*bcVals[1]+(5400.0*phiLx[0]-7875.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2182.384017536785*rdxCp2R3[1]+6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]+2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(26292.53125889555*rdxCp2[0]*rdxCp2Sq[1]+12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*volFac+((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-32735.76026305177*rdxCp2R4[1])-94154.28189944413*rdxCp2[0]*rdxCp2R3[1]-37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(2182.384017536785*rdxCp2[0]*rdxCp2R3[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(75600.0*rdxCp2R4[1]+219960.0*rdxCp2[0]*rdxCp2R3[1]+89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+28350.0*phiUy[0]*rdxCp2R4[1]+(24317.99333826703*rdxCp2[0]*phiUy[1]+36373.06695894642*rdxCp2[0]*phiLx[1]+21000.0*rdxCp2[0]*bcVals[1]+(81540.0*phiUy[0]+31500.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(11223.68923304632*rdxCp2Sq[0]*phiUy[1]+61279.95757178687*rdxCp2Sq[0]*phiLx[1]+111280.0*rdxCp2Sq[0]*bcVals[1]+(32355.0*phiUy[0]+53070.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2R3[0]*phiUy[1]+21408.14798155132*rdxCp2R3[0]*phiLx[1]+48360.0*rdxCp2R3[0]*bcVals[1]+(900.0*phiUy[0]+18540.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+692.8203230275509*rdxCp2R4[0]*phiLx[1]+1600.0*rdxCp2R4[0]*bcVals[1]+600.0*phiLx[0]*rdxCp2R4[0])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = (((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1770.0*rdxCp2[0]*rdxCp2Sq[1]+720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(31500.0*rdxCp2R3[1]+30300.0*rdxCp2[0]*rdxCp2Sq[1]+7740.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-23400.0*rdxCp2[0]*rdxCp2R3[1])-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(300.0*rdxCp2R3[0]*rdxCp2[1]-3150.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(50402.67850025433*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+28350.0*phiUy[1]*rdxCp2R4[1]+(25380.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiLx[1]+60621.7782649107*rdxCp2[0]*bcVals[1]+(20264.99444855586*phiUy[0]-45466.33369868303*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(6435.0*rdxCp2Sq[0]*phiUy[1]-25200.0*rdxCp2Sq[0]*phiLx[1]+72919.33899864974*rdxCp2Sq[0]*bcVals[1]+(9353.074360871933*phiUy[0]-21823.84017536785*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(180.0*rdxCp2R3[0]*phiUy[1]-900.0*rdxCp2R3[0]*phiLx[1]+21823.84017536785*rdxCp2R3[0]*bcVals[1]+(259.8076211353315*phiUy[0]-779.4228634059946*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+692.8203230275509*rdxCp2R4[0]*bcVals[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = (((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+5455.960043841962*rho[0]*rdxCp2R3[1]+15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiLx[3]+((-18900.0*rdxCp2R4[1])-59400.0*rdxCp2[0]*rdxCp2R3[1]-28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(5670.0*rdxCp2[0]*rdxCp2R3[1]+16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiLx[2]+((-32735.76026305177*rdxCp2R4[1])-113795.7380572752*rdxCp2[0]*rdxCp2R3[1]-64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+16367.88013152588*phiUy[0]*rdxCp2R4[1]+(6480.0*rdxCp2[0]*phiUy[1]+6300.0*rdxCp2[0]*phiLx[1]+3637.306695894642*rdxCp2[0]*bcVals[1]+(51441.90898479563*phiUy[0]+5455.960043841962*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(6480.0*rdxCp2Sq[0]*phiUy[1]+8850.0*rdxCp2Sq[0]*phiLx[1]+20438.19952931275*rdxCp2Sq[0]*bcVals[1]+(24811.62781842416*phiUy[0]+7664.324823492281*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(900.0*rdxCp2R3[0]*phiUy[1]+3000.0*rdxCp2R3[0]*phiLx[1]+7967.433714816835*rdxCp2R3[0]*bcVals[1]+(2598.076211353316*phiUy[0]+2598.076211353316*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-6235.382907247957*rdxCp2[0]*rdxCp2R3[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-8183.940065762942*rdxCp2[0]*rdxCp2R3[1])-21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]-4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(4950.0*rdxCp2[0]*rdxCp2R3[1]-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+16367.88013152588*phiUy[1]*rdxCp2R4[1]+(36476.99000740053*rdxCp2[0]*phiUy[1]-9093.266739736604*rdxCp2[0]*phiLx[1]+10500.0*rdxCp2[0]*bcVals[1]+(5400.0*phiUy[0]-7875.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(9846.708841029065*rdxCp2Sq[0]*phiUy[1]+8850.0*rdxCp2Sq[0]*bcVals[1]+5400.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2R3[0]*phiUy[1]+866.0254037844386*rdxCp2R3[0]*phiLx[1]+2600.0*rdxCp2R3[0]*bcVals[1]+(750.0*phiUy[0]+750.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(51.96152422706631*rdxCp2Sq[1]+571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+((-571.5767664977294*rdxCp2[0]*rdxCp2[1])-51.96152422706631*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*volFac+(300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(86.60254037844386*rdxCp2R3[1]+987.26896031426*rdxCp2[0]*rdxCp2Sq[1]+173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(51.96152422706631*rdxCp2[0]*rdxCp2Sq[1]+259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(200.0*rdxCp2R3[1]+2220.0*rdxCp2[0]*rdxCp2Sq[1]+100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-75.0*phiUy[0]*rdxCp2R3[1]+((-259.8076211353315*rdxCp2[0]*phiUy[1])-173.2050807568877*rdxCp2[0]*phiLx[1]-100.0*rdxCp2[0]*bcVals[1]+((-855.0*phiUy[0])-150.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-51.96152422706631*rdxCp2Sq[0]*phiUy[1])-987.26896031426*rdxCp2Sq[0]*phiLx[1]-2220.0*rdxCp2Sq[0]*bcVals[1]+((-150.0*phiUy[0])-855.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-86.60254037844386*rdxCp2R3[0]*phiLx[1]-200.0*rdxCp2R3[0]*bcVals[1]-75.0*phiLx[0]*rdxCp2R3[0]))/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]-129.9038105676658*rdxCp2Sq[0]*rho[0])*volFac+(259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(750.0*rdxCp2[0]*rdxCp2Sq[1]+150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(225.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(1991.858428704209*rdxCp2[0]*rdxCp2Sq[1]+86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-225.0*phiUy[1]*rdxCp2R3[1]+((-765.0*rdxCp2[0]*phiUy[1])+750.0*rdxCp2[0]*phiLx[1]-866.0254037844386*rdxCp2[0]*bcVals[1]+(649.5190528383289*phiLx[0]-649.5190528383289*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90.0*rdxCp2Sq[0]*phiUy[1])+150.0*rdxCp2Sq[0]*phiLx[1]-3031.088913245535*rdxCp2Sq[0]*bcVals[1]+(129.9038105676658*phiLx[0]-129.9038105676658*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]-259.8076211353315*rdxCp2R3[0]*bcVals[1]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = (((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-129.9038105676658*rho[0]*rdxCp2Sq[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiLx[3]+((-150.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(90.0*rdxCp2[0]*rdxCp2Sq[1]+765.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiLx[2]+(259.8076211353315*rdxCp2R3[1]+3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-225.0*rdxCp2[0]*phiUy[1])-150.0*rdxCp2[0]*phiLx[1]-86.60254037844386*rdxCp2[0]*bcVals[1]+(129.9038105676658*phiUy[0]-129.9038105676658*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(225.0*rdxCp2Sq[0]*phiUy[1]-750.0*rdxCp2Sq[0]*phiLx[1]-1991.858428704209*rdxCp2Sq[0]*bcVals[1]+(649.5190528383289*phiUy[0]-649.5190528383289*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+(363.7306695894642*rdxCp2[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0])*rho[2]+((-259.8076211353315*rdxCp2Sq[1])-363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])-1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(1300.0*rdxCp2[0]*rdxCp2Sq[1]+500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1299.038105676658*rdxCp2[0]*phiUy[1]+433.0127018922193*rdxCp2[0]*phiLx[1]-500.0*rdxCp2[0]*bcVals[1]+(375.0*phiLx[0]-375.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(259.8076211353315*rdxCp2Sq[0]*phiUy[1]-433.0127018922193*rdxCp2Sq[0]*phiLx[1]-1300.0*rdxCp2Sq[0]*bcVals[1]+(375.0*phiUy[0]-375.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-1454.922678357857*rdxCp2Sq[1])-8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+((-8764.17708629852*rdxCp2[0]*rdxCp2[1])-1454.922678357857*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(50400.0*rdxCp2R3[1]+314940.0*rdxCp2[0]*rdxCp2Sq[1]+63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(21823.84017536785*rdxCp2R3[1]+134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]+18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-1091.192008768392*rdxCp2[0]*rdxCp2Sq[1])-8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+18900.0*phiLy[0]*rdxCp2R3[1]+((-8105.997779422343*rdxCp2[0]*phiLy[1])+18186.53347947321*rdxCp2[0]*phiLx[1]+63000.0*rdxCp2[0]*bcVals[1]+(116685.0*phiLy[0]+15750.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-1091.192008768392*rdxCp2Sq[0]*phiLy[1])+134736.2323207829*rdxCp2Sq[0]*phiLx[1]+314940.0*rdxCp2Sq[0]*bcVals[1]+(15750.0*phiLy[0]+116685.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+21823.84017536785*rdxCp2R3[0]*phiLx[1]+50400.0*rdxCp2R3[0]*bcVals[1]+18900.0*phiLx[0]*rdxCp2R3[0])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]+3637.306695894642*rdxCp2Sq[0]*rho[0])*volFac+((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(50402.67850025433*rdxCp2[0]*rdxCp2Sq[1]+10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(23400.0*rdxCp2[0]*rdxCp2Sq[1]+3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]-18900.0*phiLy[1]*rdxCp2R3[1]+((-20745.0*rdxCp2[0]*phiLy[1])+52500.0*rdxCp2[0]*phiLx[1]-181865.3347947321*rdxCp2[0]*bcVals[1]+(20264.99444855586*phiLy[0]+45466.33369868303*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-2835.0*rdxCp2Sq[0]*phiLy[1])+85350.0*rdxCp2Sq[0]*phiLx[1]-164198.4165575295*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiLy[0]+73915.26821300182*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+12600.0*rdxCp2R3[0]*phiLx[1]-21823.84017536785*rdxCp2R3[0]*bcVals[1]+10911.92008768392*phiLx[0]*rdxCp2R3[0]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+3637.306695894642*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiLx[3]+((-21823.84017536785*rdxCp2R3[1])-164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]-181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(12600.0*rdxCp2R3[1]+85350.0*rdxCp2[0]*rdxCp2Sq[1]+52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-2835.0*rdxCp2[0]*rdxCp2Sq[1])-20745.0*rdxCp2Sq[0]*rdxCp2[1]-18900.0*rdxCp2R3[0])*phiLx[2]+10911.92008768392*phiLy[0]*rdxCp2R3[1]+((-2160.0*rdxCp2[0]*phiLy[1])+3150.0*rdxCp2[0]*phiLx[1]+10911.92008768392*rdxCp2[0]*bcVals[1]+(73915.26821300182*phiLy[0]+2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-3150.0*rdxCp2Sq[0]*phiLy[1])+23400.0*rdxCp2Sq[0]*phiLx[1]+50402.67850025433*rdxCp2Sq[0]*bcVals[1]+(45466.33369868303*phiLy[0]+20264.99444855586*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-320.4293994002423*rdxCp2[0]*rdxCp2[1])-1212.435565298214*rdxCp2Sq[0])*rho[2]+((-1212.435565298214*rdxCp2Sq[1])-320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]+1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiLx[3]+(1650.0*rdxCp2[0]*rdxCp2Sq[1]-10500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2078.460969082652*rdxCp2[0]*rdxCp2Sq[1]+3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-2727.980021920981*rdxCp2[0]*rdxCp2Sq[1])-16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]-3637.306695894642*rdxCp2R3[0])*phiLx[2]-3637.306695894642*phiLy[1]*rdxCp2R3[1]+((-16116.7327644284*rdxCp2[0]*phiLy[1])+3031.088913245535*rdxCp2[0]*phiLx[1]-10500.0*rdxCp2[0]*bcVals[1]+(1800.0*phiLy[0]+2625.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-2727.980021920981*rdxCp2Sq[0]*phiLy[1])+2078.460969082652*rdxCp2Sq[0]*phiLx[1]+1650.0*rdxCp2Sq[0]*bcVals[1]+(2625.0*phiLy[0]+1800.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-415.6921938165305*rdxCp2R3[1])-12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]-26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(2286.307065990918*rdxCp2[0]*rdxCp2Sq[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*volFac+(1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-1600.0*rdxCp2R4[1])-48360.0*rdxCp2[0]*rdxCp2R3[1]-111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-692.8203230275509*rdxCp2R4[1])-21408.14798155132*rdxCp2[0]*rdxCp2R3[1]-61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2R3[1])-11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]-24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-600.0*phiLy[0]*rdxCp2R4[1]+(1039.230484541326*rdxCp2[0]*phiLy[1]-1039.230484541326*rdxCp2[0]*phiLx[1]-3600.0*rdxCp2[0]*bcVals[1]+((-18540.0*phiLy[0])-900.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(3065.729929396912*rdxCp2Sq[0]*phiLy[1]-37360.33591926068*rdxCp2Sq[0]*phiLx[1]-89820.0*rdxCp2Sq[0]*bcVals[1]+((-53070.0*phiLy[0])-32355.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(2182.384017536785*rdxCp2R3[0]*phiLy[1]-94154.28189944413*rdxCp2R3[0]*phiLx[1]-219960.0*rdxCp2R3[0]*bcVals[1]+((-31500.0*phiLy[0])-81540.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-32735.76026305177*rdxCp2R4[0]*phiLx[1]-75600.0*rdxCp2R4[0]*bcVals[1]-28350.0*phiLx[0]*rdxCp2R4[0]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = (((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]-5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0]*rho[0])*volFac+(692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-7967.433714816835*rdxCp2[0]*rdxCp2R3[1])-20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-3000.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6480.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+600.0*phiLy[1]*rdxCp2R4[1]+(6240.0*rdxCp2[0]*phiLy[1]-3000.0*rdxCp2[0]*phiLx[1]+10392.30484541326*rdxCp2[0]*bcVals[1]+((-2598.076211353316*phiLy[0])-2598.076211353316*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(16785.0*rdxCp2Sq[0]*phiLy[1]-28650.0*rdxCp2Sq[0]*phiLx[1]+64951.90528383289*rdxCp2Sq[0]*bcVals[1]+((-7664.324823492281*phiLy[0])-24811.62781842416*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(5670.0*rdxCp2R3[0]*phiLy[1]-59400.0*rdxCp2R3[0]*phiLx[1]+113795.7380572752*rdxCp2R3[0]*bcVals[1]+((-5455.960043841962*phiLy[0])-51441.90898479563*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-18900.0*rdxCp2R4[0]*phiLx[1]+32735.76026305177*rdxCp2R4[0]*bcVals[1]-16367.88013152588*phiLx[0]*rdxCp2R4[0])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+((-240.0*rdxCp2R3[1])-7740.0*rdxCp2[0]*rdxCp2Sq[1]-30300.0*rdxCp2Sq[0]*rdxCp2[1]-31500.0*rdxCp2R3[0])*rho[2]+(720.0*rdxCp2[0]*rdxCp2Sq[1]+1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiLx[3]+((-692.8203230275509*rdxCp2R4[1])-21823.84017536785*rdxCp2[0]*rdxCp2R3[1]-72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]-60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(900.0*rdxCp2[0]*rdxCp2R3[1]+25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+52500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-180.0*rdxCp2[0]*rdxCp2R3[1])-6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25380.0*rdxCp2R3[0]*rdxCp2[1]-28350.0*rdxCp2R4[0])*phiLx[2]+(300.0*rdxCp2[0]*phiLy[1]-300.0*rdxCp2[0]*phiLx[1]-1039.230484541326*rdxCp2[0]*bcVals[1]+(779.4228634059946*phiLy[0]-259.8076211353315*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-10800.0*rdxCp2Sq[0]*phiLx[1])-24941.53162899183*rdxCp2Sq[0]*bcVals[1]+(21823.84017536785*phiLy[0]-9353.074360871933*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3150.0*rdxCp2R3[0]*phiLy[1])-23400.0*rdxCp2R3[0]*phiLx[1]-50402.67850025433*rdxCp2R3[0]*bcVals[1]+(45466.33369868303*phiLy[0]-20264.99444855586*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiLx[3]+((-2600.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(9093.266739736604*rdxCp2R3[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2R3[1])-9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]-36476.99000740053*rdxCp2R3[0]*rdxCp2[1]-16367.88013152588*rdxCp2R4[0])*phiLx[2]+((-4330.127018922193*rdxCp2[0]*phiLy[1])-866.0254037844386*rdxCp2[0]*phiLx[1]+3000.0*rdxCp2[0]*bcVals[1]+((-750.0*phiLy[0])-750.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-21823.84017536785*rdxCp2Sq[0]*phiLy[1])-6235.382907247957*rdxCp2Sq[0]*phiLx[1]+10800.0*rdxCp2Sq[0]*bcVals[1]-5400.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-8183.940065762942*rdxCp2R3[0]*phiLy[1])-6235.382907247957*rdxCp2R3[0]*phiLx[1]-4950.0*rdxCp2R3[0]*bcVals[1]+(7875.0*phiLy[0]-5400.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2182.384017536785*rdxCp2R3[1]+6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]+2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-26292.53125889555*rdxCp2[0]*rdxCp2Sq[1])-12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]-415.6921938165305*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*volFac+((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-75600.0*rdxCp2R4[1])-219960.0*rdxCp2[0]*rdxCp2R3[1]-89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-32735.76026305177*rdxCp2R4[1])-94154.28189944413*rdxCp2[0]*rdxCp2R3[1]-37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(2182.384017536785*rdxCp2[0]*rdxCp2R3[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-28350.0*phiLy[0]*rdxCp2R4[1]+((-24317.99333826703*rdxCp2[0]*phiLy[1])-36373.06695894642*rdxCp2[0]*phiLx[1]-21000.0*rdxCp2[0]*bcVals[1]+((-81540.0*phiLy[0])-31500.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-11223.68923304632*rdxCp2Sq[0]*phiLy[1])-61279.95757178687*rdxCp2Sq[0]*phiLx[1]-111280.0*rdxCp2Sq[0]*bcVals[1]+((-32355.0*phiLy[0])-53070.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2R3[0]*phiLy[1])-21408.14798155132*rdxCp2R3[0]*phiLx[1]-48360.0*rdxCp2R3[0]*bcVals[1]+((-900.0*phiLy[0])-18540.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-692.8203230275509*rdxCp2R4[0]*phiLx[1]-1600.0*rdxCp2R4[0]*bcVals[1]-600.0*phiLx[0]*rdxCp2R4[0]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1770.0*rdxCp2[0]*rdxCp2Sq[1]+720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-31500.0*rdxCp2R3[1])-30300.0*rdxCp2[0]*rdxCp2Sq[1]-7740.0*rdxCp2Sq[0]*rdxCp2[1]-240.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-50402.67850025433*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-23400.0*rdxCp2[0]*rdxCp2R3[1])-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(300.0*rdxCp2R3[0]*rdxCp2[1]-3150.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]-28350.0*phiLy[1]*rdxCp2R4[1]+((-25380.0*rdxCp2[0]*phiLy[1])+52500.0*rdxCp2[0]*phiLx[1]-60621.7782649107*rdxCp2[0]*bcVals[1]+(45466.33369868303*phiLx[0]-20264.99444855586*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-6435.0*rdxCp2Sq[0]*phiLy[1])+25200.0*rdxCp2Sq[0]*phiLx[1]-72919.33899864974*rdxCp2Sq[0]*bcVals[1]+(21823.84017536785*phiLx[0]-9353.074360871933*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-180.0*rdxCp2R3[0]*phiLy[1])+900.0*rdxCp2R3[0]*phiLx[1]-21823.84017536785*rdxCp2R3[0]*bcVals[1]+(779.4228634059946*phiLx[0]-259.8076211353315*phiLy[0])*rdxCp2R3[0])*rdxCp2[1]-692.8203230275509*rdxCp2R4[0]*bcVals[1]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = (((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-5455.960043841962*rho[0]*rdxCp2R3[1]-15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiLx[3]+(32735.76026305177*rdxCp2R4[1]+113795.7380572752*rdxCp2[0]*rdxCp2R3[1]+64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-18900.0*rdxCp2R4[1])-59400.0*rdxCp2[0]*rdxCp2R3[1]-28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(5670.0*rdxCp2[0]*rdxCp2R3[1]+16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiLx[2]-16367.88013152588*phiLy[0]*rdxCp2R4[1]+((-6480.0*rdxCp2[0]*phiLy[1])-6300.0*rdxCp2[0]*phiLx[1]-3637.306695894642*rdxCp2[0]*bcVals[1]+((-51441.90898479563*phiLy[0])-5455.960043841962*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-6480.0*rdxCp2Sq[0]*phiLy[1])-8850.0*rdxCp2Sq[0]*phiLx[1]-20438.19952931275*rdxCp2Sq[0]*bcVals[1]+((-24811.62781842416*phiLy[0])-7664.324823492281*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-900.0*rdxCp2R3[0]*phiLy[1])-3000.0*rdxCp2R3[0]*phiLx[1]-7967.433714816835*rdxCp2R3[0]*bcVals[1]+((-2598.076211353316*phiLy[0])-2598.076211353316*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-4950.0*rdxCp2[0]*rdxCp2R3[1])+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-6235.382907247957*rdxCp2[0]*rdxCp2R3[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-8183.940065762942*rdxCp2[0]*rdxCp2R3[1])-21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]-4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-16367.88013152588*phiLy[1]*rdxCp2R4[1]+((-36476.99000740053*rdxCp2[0]*phiLy[1])+9093.266739736604*rdxCp2[0]*phiLx[1]-10500.0*rdxCp2[0]*bcVals[1]+(7875.0*phiLx[0]-5400.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-9846.708841029065*rdxCp2Sq[0]*phiLy[1])-8850.0*rdxCp2Sq[0]*bcVals[1]-5400.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-519.6152422706631*rdxCp2R3[0]*phiLy[1])-866.0254037844386*rdxCp2R3[0]*phiLx[1]-2600.0*rdxCp2R3[0]*bcVals[1]+((-750.0*phiLy[0])-750.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(51.96152422706631*rdxCp2Sq[1]+571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+(571.5767664977294*rdxCp2[0]*rdxCp2[1]+51.96152422706631*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*volFac+(300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(200.0*rdxCp2R3[1]+2220.0*rdxCp2[0]*rdxCp2Sq[1]+100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(86.60254037844386*rdxCp2R3[1]+987.26896031426*rdxCp2[0]*rdxCp2Sq[1]+173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(51.96152422706631*rdxCp2[0]*rdxCp2Sq[1]+259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+75.0*phiLy[0]*rdxCp2R3[1]+(259.8076211353315*rdxCp2[0]*phiLy[1]+173.2050807568877*rdxCp2[0]*phiLx[1]+100.0*rdxCp2[0]*bcVals[1]+(855.0*phiLy[0]+150.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2Sq[0]*phiLy[1]+987.26896031426*rdxCp2Sq[0]*phiLx[1]+2220.0*rdxCp2Sq[0]*bcVals[1]+(150.0*phiLy[0]+855.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+86.60254037844386*rdxCp2R3[0]*phiLx[1]+200.0*rdxCp2R3[0]*bcVals[1]+75.0*phiLx[0]*rdxCp2R3[0])/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = (((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0]*rho[0])*volFac+(259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(1991.858428704209*rdxCp2[0]*rdxCp2Sq[1]+86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(750.0*rdxCp2[0]*rdxCp2Sq[1]+150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(225.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+225.0*phiLy[1]*rdxCp2R3[1]+(765.0*rdxCp2[0]*phiLy[1]-750.0*rdxCp2[0]*phiLx[1]+866.0254037844386*rdxCp2[0]*bcVals[1]+(649.5190528383289*phiLy[0]-649.5190528383289*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(90.0*rdxCp2Sq[0]*phiLy[1]-150.0*rdxCp2Sq[0]*phiLx[1]+3031.088913245535*rdxCp2Sq[0]*bcVals[1]+(129.9038105676658*phiLy[0]-129.9038105676658*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+259.8076211353315*rdxCp2R3[0]*bcVals[1])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = (((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+129.9038105676658*rho[0]*rdxCp2Sq[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiLx[3]+(259.8076211353315*rdxCp2R3[1]+3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-150.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(90.0*rdxCp2[0]*rdxCp2Sq[1]+765.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiLx[2]+(225.0*rdxCp2[0]*phiLy[1]+150.0*rdxCp2[0]*phiLx[1]+86.60254037844386*rdxCp2[0]*bcVals[1]+(129.9038105676658*phiLx[0]-129.9038105676658*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-225.0*rdxCp2Sq[0]*phiLy[1])+750.0*rdxCp2Sq[0]*phiLx[1]+1991.858428704209*rdxCp2Sq[0]*bcVals[1]+(649.5190528383289*phiLx[0]-649.5190528383289*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+(363.7306695894642*rdxCp2[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0])*rho[2]+(259.8076211353315*rdxCp2Sq[1]+363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(1300.0*rdxCp2[0]*rdxCp2Sq[1]+500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])-1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1299.038105676658*rdxCp2[0]*phiLy[1])-433.0127018922193*rdxCp2[0]*phiLx[1]+500.0*rdxCp2[0]*bcVals[1]+(375.0*phiLy[0]-375.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-259.8076211353315*rdxCp2Sq[0]*phiLy[1])+433.0127018922193*rdxCp2Sq[0]*phiLx[1]+1300.0*rdxCp2Sq[0]*bcVals[1]+(375.0*phiLx[0]-375.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

void MGpoissonDampedGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((6.928203230275509*rdxCp2[0]*rho[1]+12.0*rho[0]*rdxCp2[1]+100.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+5.196152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]+9.0*phiUx[0]+75.0*phiLy[0]-177.0*phiC[0]+36.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiUx[1]+(90.0*phiUx[0]-210.0*phiC[0]+240.0*bcVals[0])*rdxCp2Sq[0])*omega+18.0*phiC[0]*rdxCp2Sq[1]+177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+210.0*phiC[0]*rdxCp2Sq[0])/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[1]+34.64101615137754*rdxCp2[0]*rho[0])*omega*volFac+(((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1]-36.0*phiC[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+27.0*rdxCp2[0]*phiLy[1]-354.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUy[0]+51.96152422706631*phiUx[0]+25.98076211353316*phiLy[0]-207.8460969082653*bcVals[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiUx[1]-420.0*rdxCp2Sq[0]*phiC[1]+(103.9230484541326*phiUx[0]-207.8460969082653*bcVals[0])*rdxCp2Sq[0])*omega+36.0*phiC[1]*rdxCp2Sq[1]+354.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+420.0*rdxCp2Sq[0]*phiC[1])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = ((6.928203230275509*rdxCp2[0]*rho[3]+(80.0*rdxCp2[1]+100.0*rdxCp2[0])*rho[2])*omega*volFac+((-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[3])+((-69.28203230275508*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiUx[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-200.0*rdxCp2Sq[1])-250.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(60.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*phiUx[2]+((-200.0*rdxCp2Sq[1])-250.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-800.0*rdxCp2Sq[1])-1180.0*rdxCp2[0]*rdxCp2[1]-210.0*rdxCp2Sq[0])*phiC[2]+(173.2050807568877*phiUy[0]-173.2050807568877*phiLy[0])*rdxCp2Sq[1]+(15.0*rdxCp2[0]*phiUy[1]-15.0*rdxCp2[0]*phiLy[1]+(216.5063509461096*phiUy[0]-216.5063509461096*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(800.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0])*phiC[2])/(800.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[3] = (((160.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[0]*rho[2])*omega*volFac+(((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-400.0*rdxCp2[0]*rdxCp2[1])-120.0*rdxCp2Sq[0])*phiUx[3]+((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-1600.0*rdxCp2Sq[1])-2360.0*rdxCp2[0]*rdxCp2[1]-420.0*rdxCp2Sq[0])*phiC[3]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(346.4101615137754*rdxCp2[0]*rdxCp2[1]+103.9230484541326*rdxCp2Sq[0])*phiUx[2]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(346.4101615137754*phiUy[1]-346.4101615137754*phiLy[1])*rdxCp2Sq[1]+(77.94228634059945*rdxCp2[0]*phiUy[1]-77.94228634059945*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]-75.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0])*phiC[3])/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((41.56921938165305*rdxCp2[0]*rho[1]-36.0*rho[0]*rdxCp2[1]-120.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-27.0*phiUy[0])-27.0*phiLy[0]+54.0*phiC[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+41.56921938165305*rdxCp2[0]*phiUx[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+((-90.0*phiUy[0])-36.0*phiUx[0]-90.0*phiLy[0]+216.0*phiC[0]+24.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiUx[1]+((-60.0*phiUx[0])+60.0*phiC[0]+160.0*bcVals[0])*rdxCp2Sq[0])*omega-54.0*phiC[0]*rdxCp2Sq[1]-216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-60.0*phiC[0]*rdxCp2Sq[0]))/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[1]-34.64101615137754*rdxCp2[0]*rho[0])*omega*volFac+(((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+18.0*rdxCp2[0]*phiLy[1]-216.0*rdxCp2[0]*phiC[1]+((-25.98076211353316*phiUy[0])+51.96152422706631*phiUx[0]-25.98076211353316*phiLy[0]+69.28203230275508*bcVals[0])*rdxCp2[0])*rdxCp2[1]-60.0*rdxCp2Sq[0]*phiC[1]+69.28203230275508*bcVals[0]*rdxCp2Sq[0])*omega+54.0*phiC[1]*rdxCp2Sq[1]+216.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+60.0*rdxCp2Sq[0]*phiC[1])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((20.78460969082652*rdxCp2[0]*rho[3]+((-120.0*rdxCp2[1])-60.0*rdxCp2[0])*rho[2])*omega*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUy[3])+(138.5640646055102*rdxCp2[0]*rdxCp2[1]+34.64101615137754*rdxCp2Sq[0])*phiUx[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(300.0*rdxCp2Sq[1]+150.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-120.0*rdxCp2[0]*rdxCp2[1])-30.0*rdxCp2Sq[0])*phiUx[2]+(300.0*rdxCp2Sq[1]+150.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(1200.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0])*phiC[2]+(259.8076211353315*phiLy[0]-259.8076211353315*phiUy[0])*rdxCp2Sq[1]+(45.0*rdxCp2[0]*phiUy[1]-45.0*rdxCp2[0]*phiLy[1]+(129.9038105676658*phiLy[0]-129.9038105676658*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-1200.0*rdxCp2Sq[1])-720.0*rdxCp2[0]*rdxCp2[1]-30.0*rdxCp2Sq[0])*phiC[2]))/(1200.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0]); 
  phiC[3] = (((240.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[0]*rho[2])*omega*volFac+(((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-60.0*rdxCp2Sq[0])*phiC[3]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[2]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiLy[0]-75.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiC[3])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((6.928203230275509*rdxCp2[0]*rho[1]-12.0*rho[0]*rdxCp2[1]-100.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-9.0*phiUy[0])-9.0*phiLy[0]+18.0*phiC[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]+5.196152422706631*rdxCp2[0]*phiLy[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-36.0*rdxCp2[0]*bcVals[1]+((-75.0*phiUy[0])-75.0*phiLy[0]-9.0*phiLx[0]+177.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiLx[1]-240.0*rdxCp2Sq[0]*bcVals[1]+(210.0*phiC[0]-90.0*phiLx[0])*rdxCp2Sq[0])*omega-18.0*phiC[0]*rdxCp2Sq[1]-177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-210.0*phiC[0]*rdxCp2Sq[0]))/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[1]-34.64101615137754*rdxCp2[0]*rho[0])*omega*volFac+(((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1]-36.0*phiC[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]+27.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]-354.0*rdxCp2[0]*phiC[1]+207.8460969082653*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiUy[0])-25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiLx[1]-420.0*rdxCp2Sq[0]*phiC[1]+207.8460969082653*rdxCp2Sq[0]*bcVals[1]-103.9230484541326*phiLx[0]*rdxCp2Sq[0])*omega+36.0*phiC[1]*rdxCp2Sq[1]+354.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+420.0*rdxCp2Sq[0]*phiC[1])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((6.928203230275509*rdxCp2[0]*rho[3]+((-80.0*rdxCp2[1])-100.0*rdxCp2[0])*rho[2])*omega*volFac+((-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[3])-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-69.28203230275508*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiLx[3]+(200.0*rdxCp2Sq[1]+250.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(200.0*rdxCp2Sq[1]+250.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-60.0*rdxCp2[0]*rdxCp2[1])-90.0*rdxCp2Sq[0])*phiLx[2]+(800.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0])*phiC[2]+(173.2050807568877*phiLy[0]-173.2050807568877*phiUy[0])*rdxCp2Sq[1]+(15.0*rdxCp2[0]*phiUy[1]-15.0*rdxCp2[0]*phiLy[1]+(216.5063509461096*phiLy[0]-216.5063509461096*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-800.0*rdxCp2Sq[1])-1180.0*rdxCp2[0]*rdxCp2[1]-210.0*rdxCp2Sq[0])*phiC[2]))/(800.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[3] = (((160.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[0]*rho[2])*omega*volFac+(((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-400.0*rdxCp2Sq[1])-90.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-400.0*rdxCp2[0]*rdxCp2[1])-120.0*rdxCp2Sq[0])*phiLx[3]+((-1600.0*rdxCp2Sq[1])-2360.0*rdxCp2[0]*rdxCp2[1]-420.0*rdxCp2Sq[0])*phiC[3]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-346.4101615137754*rdxCp2[0]*rdxCp2[1])-103.9230484541326*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*phiUy[1]-346.4101615137754*phiLy[1])*rdxCp2Sq[1]+(77.94228634059945*rdxCp2[0]*phiUy[1]-77.94228634059945*rdxCp2[0]*phiLy[1]+(75.0*phiLy[0]-75.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0])*phiC[3])/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((41.56921938165305*rdxCp2[0]*rho[1]+36.0*rho[0]*rdxCp2[1]+120.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*phiUy[0]+27.0*phiLy[0]-54.0*phiC[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+41.56921938165305*rdxCp2[0]*phiLx[1]+24.0*rdxCp2[0]*bcVals[1]+(90.0*phiUy[0]+90.0*phiLy[0]+36.0*phiLx[0]-216.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiLx[1]+160.0*rdxCp2Sq[0]*bcVals[1]+(60.0*phiLx[0]-60.0*phiC[0])*rdxCp2Sq[0])*omega+54.0*phiC[0]*rdxCp2Sq[1]+216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+60.0*phiC[0]*rdxCp2Sq[0])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[1]+34.64101615137754*rdxCp2[0]*rho[0])*omega*volFac+(((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]+18.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]-216.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiUy[0]+25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]-60.0*rdxCp2Sq[0]*phiC[1]+69.28203230275508*rdxCp2Sq[0]*bcVals[1])*omega+54.0*phiC[1]*rdxCp2Sq[1]+216.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+60.0*rdxCp2Sq[0]*phiC[1])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = ((20.78460969082652*rdxCp2[0]*rho[3]+(120.0*rdxCp2[1]+60.0*rdxCp2[0])*rho[2])*omega*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUy[3])-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(138.5640646055102*rdxCp2[0]*rdxCp2[1]+34.64101615137754*rdxCp2Sq[0])*phiLx[3]+((-300.0*rdxCp2Sq[1])-150.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-300.0*rdxCp2Sq[1])-150.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(120.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0])*phiLx[2]+((-1200.0*rdxCp2Sq[1])-720.0*rdxCp2[0]*rdxCp2[1]-30.0*rdxCp2Sq[0])*phiC[2]+(259.8076211353315*phiUy[0]-259.8076211353315*phiLy[0])*rdxCp2Sq[1]+(45.0*rdxCp2[0]*phiUy[1]-45.0*rdxCp2[0]*phiLy[1]+(129.9038105676658*phiUy[0]-129.9038105676658*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1200.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0])*phiC[2])/(1200.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+30.0*rdxCp2Sq[0]); 
  phiC[3] = (((240.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[0]*rho[2])*omega*volFac+(((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-2400.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-60.0*rdxCp2Sq[0])*phiC[3]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]-346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]-75.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiC[3])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((6.928203230275509*rdxCp2[1]*rho[2]+100.0*rho[0]*rdxCp2[1]+12.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(240.0*rdxCp2Sq[1]+36.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(90.0*phiUy[0]-210.0*phiC[0])*rdxCp2Sq[1]+((-86.60254037844386*rdxCp2[0]*phiUx[1])+86.60254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUy[0]+75.0*phiUx[0]+75.0*phiLx[0]-177.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-10.39230484541326*rdxCp2Sq[0]*phiUx[1]+10.39230484541326*rdxCp2Sq[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2Sq[0])*omega+210.0*phiC[0]*rdxCp2Sq[1]+177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+18.0*phiC[0]*rdxCp2Sq[0])/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = ((6.928203230275509*rdxCp2[1]*rho[3]+(100.0*rdxCp2[1]+80.0*rdxCp2[0])*rho[1])*omega*volFac+(((-103.9230484541326*rdxCp2Sq[1])-69.28203230275508*rdxCp2[0]*rdxCp2[1])*phiUy[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUx[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLx[3]+15.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-15.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(90.0*phiUy[1]-210.0*phiC[1])*rdxCp2Sq[1]+(60.0*rdxCp2[0]*phiUy[1]-250.0*rdxCp2[0]*phiUx[1]-250.0*rdxCp2[0]*phiLx[1]-1180.0*rdxCp2[0]*phiC[1]+(216.5063509461096*phiUx[0]-216.5063509461096*phiLx[0])*rdxCp2[0])*rdxCp2[1]-200.0*rdxCp2Sq[0]*phiUx[1]-200.0*rdxCp2Sq[0]*phiLx[1]-800.0*rdxCp2Sq[0]*phiC[1]+(173.2050807568877*phiUx[0]-173.2050807568877*phiLx[0])*rdxCp2Sq[0])*omega+210.0*phiC[1]*rdxCp2Sq[1]+1180.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+800.0*rdxCp2Sq[0]*phiC[1])/(210.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+800.0*rdxCp2Sq[0]); 
  phiC[2] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[2]+34.64101615137754*rho[0]*rdxCp2[1])*omega*volFac+(((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]+((-420.0*rdxCp2Sq[1])-354.0*rdxCp2[0]*rdxCp2[1]-36.0*rdxCp2Sq[0])*phiC[2]+((-207.8460969082653*rdxCp2Sq[1])-207.8460969082653*rdxCp2[0]*rdxCp2[1])*bcVals[2]+103.9230484541326*phiUy[0]*rdxCp2Sq[1]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]+25.98076211353316*phiUx[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0])*phiC[2])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((36.0*rdxCp2[1]+160.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[1]*rho[1])*omega*volFac+(((-120.0*rdxCp2Sq[1])-400.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiUx[3]+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiLx[3]+((-420.0*rdxCp2Sq[1])-2360.0*rdxCp2[0]*rdxCp2[1]-1600.0*rdxCp2Sq[0])*phiC[3]+(77.94228634059945*rdxCp2[0]*rdxCp2[1]+346.4101615137754*rdxCp2Sq[0])*phiUx[2]+((-77.94228634059945*rdxCp2[0]*rdxCp2[1])-346.4101615137754*rdxCp2Sq[0])*phiLx[2]+103.9230484541326*phiUy[1]*rdxCp2Sq[1]+(346.4101615137754*rdxCp2[0]*phiUy[1]-86.60254037844386*rdxCp2[0]*phiUx[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiUx[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0])*phiC[3])/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((41.56921938165305*rdxCp2[1]*rho[2]-120.0*rho[0]*rdxCp2[1]-36.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiUy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(60.0*phiC[0]-60.0*phiUy[0])*rdxCp2Sq[1]+(103.9230484541326*rdxCp2[0]*phiUx[1]-103.9230484541326*rdxCp2[0]*phiLx[1]+((-36.0*phiUy[0])-90.0*phiUx[0]-90.0*phiLx[0]+216.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0]*phiUx[1]-31.17691453623978*rdxCp2Sq[0]*phiLx[1]+((-27.0*phiUx[0])-27.0*phiLx[0]+54.0*phiC[0])*rdxCp2Sq[0])*omega-60.0*phiC[0]*rdxCp2Sq[1]-216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-54.0*phiC[0]*rdxCp2Sq[0]))/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((20.78460969082652*rdxCp2[1]*rho[3]+((-60.0*rdxCp2[1])-120.0*rdxCp2[0])*rho[1])*omega*volFac+((34.64101615137754*rdxCp2Sq[1]+138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiUy[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[3]+45.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-45.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(30.0*phiC[1]-30.0*phiUy[1])*rdxCp2Sq[1]+((-120.0*rdxCp2[0]*phiUy[1])+150.0*rdxCp2[0]*phiUx[1]+150.0*rdxCp2[0]*phiLx[1]+720.0*rdxCp2[0]*phiC[1]+(129.9038105676658*phiLx[0]-129.9038105676658*phiUx[0])*rdxCp2[0])*rdxCp2[1]+300.0*rdxCp2Sq[0]*phiUx[1]+300.0*rdxCp2Sq[0]*phiLx[1]+1200.0*rdxCp2Sq[0]*phiC[1]+(259.8076211353315*phiLx[0]-259.8076211353315*phiUx[0])*rdxCp2Sq[0])*omega-30.0*phiC[1]*rdxCp2Sq[1]-720.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-1200.0*rdxCp2Sq[0]*phiC[1]))/(30.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+1200.0*rdxCp2Sq[0]); 
  phiC[2] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[2]-34.64101615137754*rho[0]*rdxCp2[1])*omega*volFac+(((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]-60.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+((-60.0*rdxCp2Sq[1])-216.0*rdxCp2[0]*rdxCp2[1]-54.0*rdxCp2Sq[0])*phiC[2]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]-25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0])*phiC[2])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((24.0*rdxCp2[1]+240.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[1]*rho[1])*omega*volFac+((-400.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+((-60.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-2400.0*rdxCp2Sq[0])*phiC[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*rdxCp2[0]*phiUy[1]+86.60254037844386*rdxCp2[0]*phiUx[1]+86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiLx[0]-75.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0])*phiC[3])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((6.928203230275509*rdxCp2[1]*rho[2]-100.0*rho[0]*rdxCp2[1]-12.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-240.0*rdxCp2Sq[1])-36.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(210.0*phiC[0]-90.0*phiLy[0])*rdxCp2Sq[1]+(86.60254037844386*rdxCp2[0]*phiUx[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+((-75.0*phiUx[0])-9.0*phiLy[0]-75.0*phiLx[0]+177.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+10.39230484541326*rdxCp2Sq[0]*phiUx[1]-10.39230484541326*rdxCp2Sq[0]*phiLx[1]+((-9.0*phiUx[0])-9.0*phiLx[0]+18.0*phiC[0])*rdxCp2Sq[0])*omega-210.0*phiC[0]*rdxCp2Sq[1]-177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-18.0*phiC[0]*rdxCp2Sq[0]))/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((6.928203230275509*rdxCp2[1]*rho[3]+((-100.0*rdxCp2[1])-80.0*rdxCp2[0])*rho[1])*omega*volFac+((-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiUx[3])+((-103.9230484541326*rdxCp2Sq[1])-69.28203230275508*rdxCp2[0]*rdxCp2[1])*phiLy[3]-17.32050807568877*rdxCp2[0]*rdxCp2[1]*phiLx[3]+15.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-15.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(210.0*phiC[1]-90.0*phiLy[1])*rdxCp2Sq[1]+(250.0*rdxCp2[0]*phiUx[1]-60.0*rdxCp2[0]*phiLy[1]+250.0*rdxCp2[0]*phiLx[1]+1180.0*rdxCp2[0]*phiC[1]+(216.5063509461096*phiLx[0]-216.5063509461096*phiUx[0])*rdxCp2[0])*rdxCp2[1]+200.0*rdxCp2Sq[0]*phiUx[1]+200.0*rdxCp2Sq[0]*phiLx[1]+800.0*rdxCp2Sq[0]*phiC[1]+(173.2050807568877*phiLx[0]-173.2050807568877*phiUx[0])*rdxCp2Sq[0])*omega-210.0*phiC[1]*rdxCp2Sq[1]-1180.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-800.0*rdxCp2Sq[0]*phiC[1]))/(210.0*rdxCp2Sq[1]+1180.0*rdxCp2[0]*rdxCp2[1]+800.0*rdxCp2Sq[0]); 
  phiC[2] = (((36.0*rdxCp2[1]+24.0*rdxCp2[0])*rho[2]-34.64101615137754*rho[0]*rdxCp2[1])*omega*volFac+(((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+(207.8460969082653*rdxCp2Sq[1]+207.8460969082653*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]+((-420.0*rdxCp2Sq[1])-354.0*rdxCp2[0]*rdxCp2[1]-36.0*rdxCp2Sq[0])*phiC[2]-103.9230484541326*phiLy[0]*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+((-25.98076211353316*phiUx[0])-51.96152422706631*phiLy[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0])*phiC[2])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((36.0*rdxCp2[1]+160.0*rdxCp2[0])*rho[3]-34.64101615137754*rdxCp2[1]*rho[1])*omega*volFac+(((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiUx[3]+((-120.0*rdxCp2Sq[1])-400.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-90.0*rdxCp2[0]*rdxCp2[1])-400.0*rdxCp2Sq[0])*phiLx[3]+((-420.0*rdxCp2Sq[1])-2360.0*rdxCp2[0]*rdxCp2[1]-1600.0*rdxCp2Sq[0])*phiC[3]+(77.94228634059945*rdxCp2[0]*rdxCp2[1]+346.4101615137754*rdxCp2Sq[0])*phiUx[2]+((-77.94228634059945*rdxCp2[0]*rdxCp2[1])-346.4101615137754*rdxCp2Sq[0])*phiLx[2]-103.9230484541326*phiLy[1]*rdxCp2Sq[1]+(86.60254037844386*rdxCp2[0]*phiUx[1]-346.4101615137754*rdxCp2[0]*phiLy[1]+86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiLx[0]-75.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0])*phiC[3])/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((41.56921938165305*rdxCp2[1]*rho[2]+120.0*rho[0]*rdxCp2[1]+36.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiLy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(60.0*phiLy[0]-60.0*phiC[0])*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2[0]*phiUx[1])+103.9230484541326*rdxCp2[0]*phiLx[1]+(90.0*phiUx[0]+36.0*phiLy[0]+90.0*phiLx[0]-216.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-31.17691453623978*rdxCp2Sq[0]*phiUx[1]+31.17691453623978*rdxCp2Sq[0]*phiLx[1]+(27.0*phiUx[0]+27.0*phiLx[0]-54.0*phiC[0])*rdxCp2Sq[0])*omega+60.0*phiC[0]*rdxCp2Sq[1]+216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+54.0*phiC[0]*rdxCp2Sq[0])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = ((20.78460969082652*rdxCp2[1]*rho[3]+(60.0*rdxCp2[1]+120.0*rdxCp2[0])*rho[1])*omega*volFac+((-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[3])+(34.64101615137754*rdxCp2Sq[1]+138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiLy[3]-51.96152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[3]+45.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-45.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(30.0*phiLy[1]-30.0*phiC[1])*rdxCp2Sq[1]+((-150.0*rdxCp2[0]*phiUx[1])+120.0*rdxCp2[0]*phiLy[1]-150.0*rdxCp2[0]*phiLx[1]-720.0*rdxCp2[0]*phiC[1]+(129.9038105676658*phiUx[0]-129.9038105676658*phiLx[0])*rdxCp2[0])*rdxCp2[1]-300.0*rdxCp2Sq[0]*phiUx[1]-300.0*rdxCp2Sq[0]*phiLx[1]-1200.0*rdxCp2Sq[0]*phiC[1]+(259.8076211353315*phiUx[0]-259.8076211353315*phiLx[0])*rdxCp2Sq[0])*omega+30.0*phiC[1]*rdxCp2Sq[1]+720.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+1200.0*rdxCp2Sq[0]*phiC[1])/(30.0*rdxCp2Sq[1]+720.0*rdxCp2[0]*rdxCp2[1]+1200.0*rdxCp2Sq[0]); 
  phiC[2] = (((24.0*rdxCp2[1]+36.0*rdxCp2[0])*rho[2]+34.64101615137754*rho[0]*rdxCp2[1])*omega*volFac+(((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]-60.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+((-60.0*rdxCp2Sq[1])-216.0*rdxCp2[0]*rdxCp2[1]-54.0*rdxCp2Sq[0])*phiC[2]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(25.98076211353316*phiUx[0]-51.96152422706631*phiLy[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0])*phiC[2])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((24.0*rdxCp2[1]+240.0*rdxCp2[0])*rho[3]+34.64101615137754*rdxCp2[1]*rho[1])*omega*volFac+(((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+((-60.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-2400.0*rdxCp2Sq[0])*phiC[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+((-86.60254037844386*rdxCp2[0]*phiUx[1])-346.4101615137754*rdxCp2[0]*phiLy[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiUx[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0])*phiC[3])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(1454.922678357857*rdxCp2Sq[1]+8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+(8764.17708629852*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-21823.84017536785*rdxCp2R3[1])-134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]-18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(1091.192008768392*rdxCp2[0]*rdxCp2Sq[1]+8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(50400.0*rdxCp2R3[1]+314940.0*rdxCp2[0]*rdxCp2Sq[1]+63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(18900.0*phiUy[0]-44100.0*phiC[0])*rdxCp2R3[1]+(8105.997779422343*rdxCp2[0]*phiUy[1]-18186.53347947321*rdxCp2[0]*phiUx[1]+(116685.0*phiUy[0]+15750.0*phiUx[0]-321405.0*phiC[0]+63000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1091.192008768392*rdxCp2Sq[0]*phiUy[1]-134736.2323207829*rdxCp2Sq[0]*phiUx[1]+(15750.0*phiUy[0]+116685.0*phiUx[0]-321405.0*phiC[0]+314940.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-21823.84017536785*rdxCp2R3[0]*phiUx[1]+(18900.0*phiUx[0]-44100.0*phiC[0]+50400.0*bcVals[0])*rdxCp2R3[0])*omega+44100.0*phiC[0]*rdxCp2R3[1]+321405.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+321405.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+44100.0*phiC[0]*rdxCp2R3[0])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = (((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]+3637.306695894642*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-23400.0*rdxCp2[0]*rdxCp2Sq[1])-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(3150.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(50402.67850025433*rdxCp2[0]*rdxCp2Sq[1]+10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(18900.0*phiUy[1]-44100.0*phiC[1])*rdxCp2R3[1]+(20745.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiUx[1]-321405.0*rdxCp2[0]*phiC[1]+(20264.99444855586*phiUy[0]+45466.33369868303*phiUx[0]-181865.3347947321*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2835.0*rdxCp2Sq[0]*phiUy[1]-85350.0*rdxCp2Sq[0]*phiUx[1]-321405.0*rdxCp2Sq[0]*phiC[1]+(2727.980021920981*phiUy[0]+73915.26821300182*phiUx[0]-164198.4165575295*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-12600.0*rdxCp2R3[0]*phiUx[1]-44100.0*rdxCp2R3[0]*phiC[1]+(10911.92008768392*phiUx[0]-21823.84017536785*bcVals[0])*rdxCp2R3[0])*omega+44100.0*phiC[1]*rdxCp2R3[1]+321405.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+44100.0*rdxCp2R3[0]*phiC[1])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = (((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+3637.306695894642*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiUx[3]+((-12600.0*rdxCp2R3[1])-85350.0*rdxCp2[0]*rdxCp2Sq[1]-52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(2835.0*rdxCp2[0]*rdxCp2Sq[1]+20745.0*rdxCp2Sq[0]*rdxCp2[1]+18900.0*rdxCp2R3[0])*phiUx[2]+((-44100.0*rdxCp2R3[1])-321405.0*rdxCp2[0]*rdxCp2Sq[1]-321405.0*rdxCp2Sq[0]*rdxCp2[1]-44100.0*rdxCp2R3[0])*phiC[2]+((-21823.84017536785*rdxCp2R3[1])-164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]-181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+10911.92008768392*phiUy[0]*rdxCp2R3[1]+(2160.0*rdxCp2[0]*phiUy[1]-3150.0*rdxCp2[0]*phiUx[1]+(73915.26821300182*phiUy[0]+2727.980021920981*phiUx[0]+10911.92008768392*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(3150.0*rdxCp2Sq[0]*phiUy[1]-23400.0*rdxCp2Sq[0]*phiUx[1]+(45466.33369868303*phiUy[0]+20264.99444855586*phiUx[0]+50402.67850025433*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0])*phiC[2])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(320.4293994002423*rdxCp2[0]*rdxCp2[1]+1212.435565298214*rdxCp2Sq[0])*rho[2]+(1212.435565298214*rdxCp2Sq[1]+320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]+1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiUx[3]+((-14700.0*rdxCp2R3[1])-107135.0*rdxCp2[0]*rdxCp2Sq[1]-107135.0*rdxCp2Sq[0]*rdxCp2[1]-14700.0*rdxCp2R3[0])*phiC[3]+((-2078.460969082652*rdxCp2[0]*rdxCp2Sq[1])-3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(2727.980021920981*rdxCp2[0]*rdxCp2Sq[1]+16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]+3637.306695894642*rdxCp2R3[0])*phiUx[2]+(1650.0*rdxCp2[0]*rdxCp2Sq[1]-10500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+3637.306695894642*phiUy[1]*rdxCp2R3[1]+(16116.7327644284*rdxCp2[0]*phiUy[1]-3031.088913245535*rdxCp2[0]*phiUx[1]+(1800.0*phiUy[0]+2625.0*phiUx[0]-10500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2727.980021920981*rdxCp2Sq[0]*phiUy[1]-2078.460969082652*rdxCp2Sq[0]*phiUx[1]+(2625.0*phiUy[0]+1800.0*phiUx[0]+1650.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0])*phiC[3])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(415.6921938165305*rdxCp2R3[1]+12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]+26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-2286.307065990918*rdxCp2[0]*rdxCp2Sq[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]-2182.384017536785*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+((1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(692.8203230275509*rdxCp2R4[1]+21408.14798155132*rdxCp2[0]*rdxCp2R3[1]+61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(311.7691453623978*rdxCp2[0]*rdxCp2R3[1]+11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]+24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(1600.0*rdxCp2R4[1]+48360.0*rdxCp2[0]*rdxCp2R3[1]+111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(600.0*phiC[0]-600.0*phiUy[0])*rdxCp2R4[1]+((-1039.230484541326*rdxCp2[0]*phiUy[1])+1039.230484541326*rdxCp2[0]*phiUx[1]+((-18540.0*phiUy[0])-900.0*phiUx[0]+21240.0*phiC[0]-3600.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-3065.729929396912*rdxCp2Sq[0]*phiUy[1])+37360.33591926068*rdxCp2Sq[0]*phiUx[1]+((-53070.0*phiUy[0])-32355.0*phiUx[0]+130335.0*phiC[0]-89820.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2182.384017536785*rdxCp2R3[0]*phiUy[1])+94154.28189944413*rdxCp2R3[0]*phiUx[1]+((-31500.0*phiUy[0])-81540.0*phiUx[0]+223020.0*phiC[0]-219960.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+32735.76026305177*rdxCp2R4[0]*phiUx[1]+((-28350.0*phiUx[0])+66150.0*phiC[0]-75600.0*bcVals[0])*rdxCp2R4[0])*omega-600.0*phiC[0]*rdxCp2R4[1]-21240.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-66150.0*phiC[0]*rdxCp2R4[0]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]-5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0]*rho[0])*omega*volFac+((692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(3000.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(900.0*rdxCp2[0]*rdxCp2R3[1]+6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6480.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(7967.433714816835*rdxCp2[0]*rdxCp2R3[1]+20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]+3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(600.0*phiC[1]-600.0*phiUy[1])*rdxCp2R4[1]+((-6240.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiUx[1]+21240.0*rdxCp2[0]*phiC[1]+((-2598.076211353316*phiUy[0])-2598.076211353316*phiUx[0]+10392.30484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-16785.0*rdxCp2Sq[0]*phiUy[1])+28650.0*rdxCp2Sq[0]*phiUx[1]+130335.0*rdxCp2Sq[0]*phiC[1]+((-7664.324823492281*phiUy[0])-24811.62781842416*phiUx[0]+64951.90528383289*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-5670.0*rdxCp2R3[0]*phiUy[1])+59400.0*rdxCp2R3[0]*phiUx[1]+223020.0*rdxCp2R3[0]*phiC[1]+((-5455.960043841962*phiUy[0])-51441.90898479563*phiUx[0]+113795.7380572752*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+18900.0*rdxCp2R4[0]*phiUx[1]+66150.0*rdxCp2R4[0]*phiC[1]+(32735.76026305177*bcVals[0]-16367.88013152588*phiUx[0])*rdxCp2R4[0])*omega-600.0*phiC[1]*rdxCp2R4[1]-21240.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-66150.0*rdxCp2R4[0]*phiC[1]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = (((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+(240.0*rdxCp2R3[1]+7740.0*rdxCp2[0]*rdxCp2Sq[1]+30300.0*rdxCp2Sq[0]*rdxCp2[1]+31500.0*rdxCp2R3[0])*rho[2]+((-720.0*rdxCp2[0]*rdxCp2Sq[1])-1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiUx[3]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-52500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(180.0*rdxCp2[0]*rdxCp2R3[1]+6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25380.0*rdxCp2R3[0]*rdxCp2[1]+28350.0*rdxCp2R4[0])*phiUx[2]+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[2]+(692.8203230275509*rdxCp2R4[1]+21823.84017536785*rdxCp2[0]*rdxCp2R3[1]+72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]+60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-300.0*rdxCp2[0]*phiUy[1])+300.0*rdxCp2[0]*phiUx[1]+(779.4228634059946*phiUy[0]-259.8076211353315*phiUx[0]-1039.230484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(10800.0*rdxCp2Sq[0]*phiUx[1]+(21823.84017536785*phiUy[0]-9353.074360871933*phiUx[0]-24941.53162899183*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3150.0*rdxCp2R3[0]*phiUy[1]+23400.0*rdxCp2R3[0]*phiUx[1]+(45466.33369868303*phiUy[0]-20264.99444855586*phiUx[0]-50402.67850025433*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[2])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiUx[3]+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[3]+(866.0254037844386*rdxCp2[0]*rdxCp2R3[1]-9093.266739736604*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(519.6152422706631*rdxCp2[0]*rdxCp2R3[1]+9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]+36476.99000740053*rdxCp2R3[0]*rdxCp2[1]+16367.88013152588*rdxCp2R4[0])*phiUx[2]+(2600.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(4330.127018922193*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiUx[1]+((-750.0*phiUy[0])-750.0*phiUx[0]+3000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]+6235.382907247957*rdxCp2Sq[0]*phiUx[1]+(10800.0*bcVals[0]-5400.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(8183.940065762942*rdxCp2R3[0]*phiUy[1]+6235.382907247957*rdxCp2R3[0]*phiUx[1]+(7875.0*phiUy[0]-5400.0*phiUx[0]-4950.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[3])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2182.384017536785*rdxCp2R3[1])-6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]-2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(26292.53125889555*rdxCp2[0]*rdxCp2Sq[1]+12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(32735.76026305177*rdxCp2R4[1]+94154.28189944413*rdxCp2[0]*rdxCp2R3[1]+37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-2182.384017536785*rdxCp2[0]*rdxCp2R3[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-75600.0*rdxCp2R4[1])-219960.0*rdxCp2[0]*rdxCp2R3[1]-89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(66150.0*phiC[0]-28350.0*phiUy[0])*rdxCp2R4[1]+(24317.99333826703*rdxCp2[0]*phiUy[1]+36373.06695894642*rdxCp2[0]*phiUx[1]+((-81540.0*phiUy[0])-31500.0*phiUx[0]+223020.0*phiC[0]+21000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(11223.68923304632*rdxCp2Sq[0]*phiUy[1]+61279.95757178687*rdxCp2Sq[0]*phiUx[1]+((-32355.0*phiUy[0])-53070.0*phiUx[0]+130335.0*phiC[0]+111280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2R3[0]*phiUy[1]+21408.14798155132*rdxCp2R3[0]*phiUx[1]+((-900.0*phiUy[0])-18540.0*phiUx[0]+21240.0*phiC[0]+48360.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+692.8203230275509*rdxCp2R4[0]*phiUx[1]+((-600.0*phiUx[0])+600.0*phiC[0]+1600.0*bcVals[0])*rdxCp2R4[0])*omega-66150.0*phiC[0]*rdxCp2R4[1]-223020.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-600.0*phiC[0]*rdxCp2R4[0]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = (((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1770.0*rdxCp2[0]*rdxCp2Sq[1])-720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(31500.0*rdxCp2R3[1]+30300.0*rdxCp2[0]*rdxCp2Sq[1]+7740.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+(((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(23400.0*rdxCp2[0]*rdxCp2R3[1]+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3150.0*rdxCp2[0]*rdxCp2R3[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-50402.67850025433*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(28350.0*phiUy[1]-66150.0*phiC[1])*rdxCp2R4[1]+(25380.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiUx[1]-223020.0*rdxCp2[0]*phiC[1]+((-20264.99444855586*phiUy[0])+45466.33369868303*phiUx[0]+60621.7782649107*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(6435.0*rdxCp2Sq[0]*phiUy[1]-25200.0*rdxCp2Sq[0]*phiUx[1]-130335.0*rdxCp2Sq[0]*phiC[1]+((-9353.074360871933*phiUy[0])+21823.84017536785*phiUx[0]+72919.33899864974*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(180.0*rdxCp2R3[0]*phiUy[1]-900.0*rdxCp2R3[0]*phiUx[1]-21240.0*rdxCp2R3[0]*phiC[1]+((-259.8076211353315*phiUy[0])+779.4228634059946*phiUx[0]+21823.84017536785*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-600.0*rdxCp2R4[0]*phiC[1]+692.8203230275509*bcVals[0]*rdxCp2R4[0])*omega+66150.0*phiC[1]*rdxCp2R4[1]+223020.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+600.0*rdxCp2R4[0]*phiC[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-5455.960043841962*rho[0]*rdxCp2R3[1]-15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiUx[3]+(18900.0*rdxCp2R4[1]+59400.0*rdxCp2[0]*rdxCp2R3[1]+28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-5670.0*rdxCp2[0]*rdxCp2R3[1])-16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiUx[2]+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[2]+(32735.76026305177*rdxCp2R4[1]+113795.7380572752*rdxCp2[0]*rdxCp2R3[1]+64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-16367.88013152588*phiUy[0]*rdxCp2R4[1]+(6480.0*rdxCp2[0]*phiUy[1]+6300.0*rdxCp2[0]*phiUx[1]+((-51441.90898479563*phiUy[0])-5455.960043841962*phiUx[0]+3637.306695894642*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(6480.0*rdxCp2Sq[0]*phiUy[1]+8850.0*rdxCp2Sq[0]*phiUx[1]+((-24811.62781842416*phiUy[0])-7664.324823492281*phiUx[0]+20438.19952931275*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(900.0*rdxCp2R3[0]*phiUy[1]+3000.0*rdxCp2R3[0]*phiUx[1]+((-2598.076211353316*phiUy[0])-2598.076211353316*phiUx[0]+7967.433714816835*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[2]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[3]+(6235.382907247957*rdxCp2[0]*rdxCp2R3[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(8183.940065762942*rdxCp2[0]*rdxCp2R3[1]+21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]+4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-4950.0*rdxCp2[0]*rdxCp2R3[1])+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+16367.88013152588*phiUy[1]*rdxCp2R4[1]+(36476.99000740053*rdxCp2[0]*phiUy[1]-9093.266739736604*rdxCp2[0]*phiUx[1]+((-5400.0*phiUy[0])+7875.0*phiUx[0]+10500.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(9846.708841029065*rdxCp2Sq[0]*phiUy[1]+(8850.0*bcVals[0]-5400.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2R3[0]*phiUy[1]+866.0254037844386*rdxCp2R3[0]*phiUx[1]+((-750.0*phiUy[0])-750.0*phiUx[0]+2600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[3])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-51.96152422706631*rdxCp2Sq[1])-571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+((-571.5767664977294*rdxCp2[0]*rdxCp2[1])-51.96152422706631*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-86.60254037844386*rdxCp2R3[1])-987.26896031426*rdxCp2[0]*rdxCp2Sq[1]-173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2Sq[1])-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-200.0*rdxCp2R3[1])-2220.0*rdxCp2[0]*rdxCp2Sq[1]-100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(75.0*phiUy[0]-75.0*phiC[0])*rdxCp2R3[1]+((-259.8076211353315*rdxCp2[0]*phiUy[1])-173.2050807568877*rdxCp2[0]*phiUx[1]+(855.0*phiUy[0]+150.0*phiUx[0]-1005.0*phiC[0]-100.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-51.96152422706631*rdxCp2Sq[0]*phiUy[1])-987.26896031426*rdxCp2Sq[0]*phiUx[1]+(150.0*phiUy[0]+855.0*phiUx[0]-1005.0*phiC[0]-2220.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-86.60254037844386*rdxCp2R3[0]*phiUx[1]+(75.0*phiUx[0]-75.0*phiC[0]-200.0*bcVals[0])*rdxCp2R3[0])*omega+75.0*phiC[0]*rdxCp2R3[1]+1005.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+1005.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+75.0*phiC[0]*rdxCp2R3[0])/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0]*rho[0])*omega*volFac+((259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(225.0*rdxCp2[0]*rdxCp2Sq[1]-225.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1991.858428704209*rdxCp2[0]*rdxCp2Sq[1])-86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(225.0*phiC[1]-225.0*phiUy[1])*rdxCp2R3[1]+((-765.0*rdxCp2[0]*phiUy[1])+750.0*rdxCp2[0]*phiUx[1]+3015.0*rdxCp2[0]*phiC[1]+(649.5190528383289*phiUy[0]-649.5190528383289*phiUx[0]-866.0254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90.0*rdxCp2Sq[0]*phiUy[1])+150.0*rdxCp2Sq[0]*phiUx[1]+3015.0*rdxCp2Sq[0]*phiC[1]+(129.9038105676658*phiUy[0]-129.9038105676658*phiUx[0]-3031.088913245535*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+225.0*rdxCp2R3[0]*phiC[1]-259.8076211353315*bcVals[0]*rdxCp2R3[0])*omega-225.0*phiC[1]*rdxCp2R3[1]-3015.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-3015.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-225.0*rdxCp2R3[0]*phiC[1]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+129.9038105676658*rho[0]*rdxCp2Sq[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiUx[3]+(150.0*rdxCp2[0]*rdxCp2Sq[1]+750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-90.0*rdxCp2[0]*rdxCp2Sq[1])-765.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiUx[2]+(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiC[2]+((-259.8076211353315*rdxCp2R3[1])-3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-225.0*rdxCp2[0]*phiUy[1])-150.0*rdxCp2[0]*phiUx[1]+((-129.9038105676658*phiUy[0])+129.9038105676658*phiUx[0]-86.60254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(225.0*rdxCp2Sq[0]*phiUy[1]-750.0*rdxCp2Sq[0]*phiUx[1]+((-649.5190528383289*phiUy[0])+649.5190528383289*phiUx[0]-1991.858428704209*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-225.0*rdxCp2R3[1])-3015.0*rdxCp2[0]*rdxCp2Sq[1]-3015.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiC[2]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+((-363.7306695894642*rdxCp2[0]*rdxCp2[1])-259.8076211353315*rdxCp2Sq[0])*rho[2]+((-259.8076211353315*rdxCp2Sq[1])-363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-450.0*rdxCp2R3[1])-6030.0*rdxCp2[0]*rdxCp2Sq[1]-6030.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2R3[0])*phiC[3]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]+1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-1300.0*rdxCp2[0]*rdxCp2Sq[1])-500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1299.038105676658*rdxCp2[0]*phiUy[1]+433.0127018922193*rdxCp2[0]*phiUx[1]+(375.0*phiUy[0]-375.0*phiUx[0]-500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(259.8076211353315*rdxCp2Sq[0]*phiUy[1]-433.0127018922193*rdxCp2Sq[0]*phiUx[1]+((-375.0*phiUy[0])+375.0*phiUx[0]-1300.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0])*phiC[3])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(1454.922678357857*rdxCp2Sq[1]+8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+((-8764.17708629852*rdxCp2[0]*rdxCp2[1])-1454.922678357857*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-50400.0*rdxCp2R3[1])-314940.0*rdxCp2[0]*rdxCp2Sq[1]-63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(1091.192008768392*rdxCp2[0]*rdxCp2Sq[1]+8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-21823.84017536785*rdxCp2R3[1])-134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]-18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(44100.0*phiC[0]-18900.0*phiLy[0])*rdxCp2R3[1]+(18186.53347947321*rdxCp2[0]*phiUx[1]-8105.997779422343*rdxCp2[0]*phiLy[1]+((-15750.0*phiUx[0])-116685.0*phiLy[0]+321405.0*phiC[0]-63000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(134736.2323207829*rdxCp2Sq[0]*phiUx[1]-1091.192008768392*rdxCp2Sq[0]*phiLy[1]+((-116685.0*phiUx[0])-15750.0*phiLy[0]+321405.0*phiC[0]-314940.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+21823.84017536785*rdxCp2R3[0]*phiUx[1]+((-18900.0*phiUx[0])+44100.0*phiC[0]-50400.0*bcVals[0])*rdxCp2R3[0])*omega-44100.0*phiC[0]*rdxCp2R3[1]-321405.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-321405.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-44100.0*phiC[0]*rdxCp2R3[0]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]-3637.306695894642*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-50402.67850025433*rdxCp2[0]*rdxCp2Sq[1])-10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(3150.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-23400.0*rdxCp2[0]*rdxCp2Sq[1])-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(44100.0*phiC[1]-18900.0*phiLy[1])*rdxCp2R3[1]+(52500.0*rdxCp2[0]*phiUx[1]-20745.0*rdxCp2[0]*phiLy[1]+321405.0*rdxCp2[0]*phiC[1]+((-45466.33369868303*phiUx[0])-20264.99444855586*phiLy[0]+181865.3347947321*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(85350.0*rdxCp2Sq[0]*phiUx[1]-2835.0*rdxCp2Sq[0]*phiLy[1]+321405.0*rdxCp2Sq[0]*phiC[1]+((-73915.26821300182*phiUx[0])-2727.980021920981*phiLy[0]+164198.4165575295*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+12600.0*rdxCp2R3[0]*phiUx[1]+44100.0*rdxCp2R3[0]*phiC[1]+(21823.84017536785*bcVals[0]-10911.92008768392*phiUx[0])*rdxCp2R3[0])*omega-44100.0*phiC[1]*rdxCp2R3[1]-321405.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-321405.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-44100.0*rdxCp2R3[0]*phiC[1]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = (((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-3637.306695894642*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiUx[3]+((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(21823.84017536785*rdxCp2R3[1]+164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]+181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2835.0*rdxCp2[0]*rdxCp2Sq[1]+20745.0*rdxCp2Sq[0]*rdxCp2[1]+18900.0*rdxCp2R3[0])*phiUx[2]+((-12600.0*rdxCp2R3[1])-85350.0*rdxCp2[0]*rdxCp2Sq[1]-52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-44100.0*rdxCp2R3[1])-321405.0*rdxCp2[0]*rdxCp2Sq[1]-321405.0*rdxCp2Sq[0]*rdxCp2[1]-44100.0*rdxCp2R3[0])*phiC[2]-10911.92008768392*phiLy[0]*rdxCp2R3[1]+(3150.0*rdxCp2[0]*phiUx[1]-2160.0*rdxCp2[0]*phiLy[1]+((-2727.980021920981*phiUx[0])-73915.26821300182*phiLy[0]-10911.92008768392*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(23400.0*rdxCp2Sq[0]*phiUx[1]-3150.0*rdxCp2Sq[0]*phiLy[1]+((-20264.99444855586*phiUx[0])-45466.33369868303*phiLy[0]-50402.67850025433*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0])*phiC[2])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(320.4293994002423*rdxCp2[0]*rdxCp2[1]+1212.435565298214*rdxCp2Sq[0])*rho[2]+((-1212.435565298214*rdxCp2Sq[1])-320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]-1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiUx[3]+((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-14700.0*rdxCp2R3[1])-107135.0*rdxCp2[0]*rdxCp2Sq[1]-107135.0*rdxCp2Sq[0]*rdxCp2[1]-14700.0*rdxCp2R3[0])*phiC[3]+(10500.0*rdxCp2Sq[0]*rdxCp2[1]-1650.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(2727.980021920981*rdxCp2[0]*rdxCp2Sq[1]+16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]+3637.306695894642*rdxCp2R3[0])*phiUx[2]+((-2078.460969082652*rdxCp2[0]*rdxCp2Sq[1])-3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-3637.306695894642*phiLy[1]*rdxCp2R3[1]+(3031.088913245535*rdxCp2[0]*phiUx[1]-16116.7327644284*rdxCp2[0]*phiLy[1]+((-2625.0*phiUx[0])-1800.0*phiLy[0]+10500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2078.460969082652*rdxCp2Sq[0]*phiUx[1]-2727.980021920981*rdxCp2Sq[0]*phiLy[1]+((-1800.0*phiUx[0])-2625.0*phiLy[0]-1650.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0])*phiC[3])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(415.6921938165305*rdxCp2R3[1]+12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]+26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(2286.307065990918*rdxCp2[0]*rdxCp2Sq[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(1600.0*rdxCp2R4[1]+48360.0*rdxCp2[0]*rdxCp2R3[1]+111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(311.7691453623978*rdxCp2[0]*rdxCp2R3[1]+11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]+24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(692.8203230275509*rdxCp2R4[1]+21408.14798155132*rdxCp2[0]*rdxCp2R3[1]+61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(600.0*phiLy[0]-600.0*phiC[0])*rdxCp2R4[1]+((-1039.230484541326*rdxCp2[0]*phiUx[1])+1039.230484541326*rdxCp2[0]*phiLy[1]+(900.0*phiUx[0]+18540.0*phiLy[0]-21240.0*phiC[0]+3600.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-37360.33591926068*rdxCp2Sq[0]*phiUx[1])+3065.729929396912*rdxCp2Sq[0]*phiLy[1]+(32355.0*phiUx[0]+53070.0*phiLy[0]-130335.0*phiC[0]+89820.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-94154.28189944413*rdxCp2R3[0]*phiUx[1])+2182.384017536785*rdxCp2R3[0]*phiLy[1]+(81540.0*phiUx[0]+31500.0*phiLy[0]-223020.0*phiC[0]+219960.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-32735.76026305177*rdxCp2R4[0]*phiUx[1]+(28350.0*phiUx[0]-66150.0*phiC[0]+75600.0*bcVals[0])*rdxCp2R4[0])*omega+600.0*phiC[0]*rdxCp2R4[1]+21240.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+66150.0*phiC[0]*rdxCp2R4[0])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = (((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]+5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0]*rho[0])*omega*volFac+(((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(7967.433714816835*rdxCp2[0]*rdxCp2R3[1]+20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]+3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(900.0*rdxCp2[0]*rdxCp2R3[1]+6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6480.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(3000.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(600.0*phiLy[1]-600.0*phiC[1])*rdxCp2R4[1]+((-3000.0*rdxCp2[0]*phiUx[1])+6240.0*rdxCp2[0]*phiLy[1]-21240.0*rdxCp2[0]*phiC[1]+(2598.076211353316*phiUx[0]+2598.076211353316*phiLy[0]-10392.30484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-28650.0*rdxCp2Sq[0]*phiUx[1])+16785.0*rdxCp2Sq[0]*phiLy[1]-130335.0*rdxCp2Sq[0]*phiC[1]+(24811.62781842416*phiUx[0]+7664.324823492281*phiLy[0]-64951.90528383289*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-59400.0*rdxCp2R3[0]*phiUx[1])+5670.0*rdxCp2R3[0]*phiLy[1]-223020.0*rdxCp2R3[0]*phiC[1]+(51441.90898479563*phiUx[0]+5455.960043841962*phiLy[0]-113795.7380572752*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-18900.0*rdxCp2R4[0]*phiUx[1]-66150.0*rdxCp2R4[0]*phiC[1]+(16367.88013152588*phiUx[0]-32735.76026305177*bcVals[0])*rdxCp2R4[0])*omega+600.0*phiC[1]*rdxCp2R4[1]+21240.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+66150.0*rdxCp2R4[0]*phiC[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = (((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+(240.0*rdxCp2R3[1]+7740.0*rdxCp2[0]*rdxCp2Sq[1]+30300.0*rdxCp2Sq[0]*rdxCp2[1]+31500.0*rdxCp2R3[0])*rho[2]+(720.0*rdxCp2[0]*rdxCp2Sq[1]+1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiUx[3]+(346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(692.8203230275509*rdxCp2R4[1]+21823.84017536785*rdxCp2[0]*rdxCp2R3[1]+72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]+60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(180.0*rdxCp2[0]*rdxCp2R3[1]+6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25380.0*rdxCp2R3[0]*rdxCp2[1]+28350.0*rdxCp2R4[0])*phiUx[2]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-52500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[2]+((-300.0*rdxCp2[0]*phiUx[1])+300.0*rdxCp2[0]*phiLy[1]+(259.8076211353315*phiUx[0]-779.4228634059946*phiLy[0]+1039.230484541326*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((9353.074360871933*phiUx[0]-21823.84017536785*phiLy[0]+24941.53162899183*bcVals[0])*rdxCp2Sq[0]-10800.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-23400.0*rdxCp2R3[0]*phiUx[1])-3150.0*rdxCp2R3[0]*phiLy[1]+(20264.99444855586*phiUx[0]-45466.33369868303*phiLy[0]+50402.67850025433*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[2])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiUx[3]+((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[3]+(2600.0*rdxCp2[0]*rdxCp2R3[1]+8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]+10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(519.6152422706631*rdxCp2[0]*rdxCp2R3[1]+9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]+36476.99000740053*rdxCp2R3[0]*rdxCp2[1]+16367.88013152588*rdxCp2R4[0])*phiUx[2]+(866.0254037844386*rdxCp2[0]*rdxCp2R3[1]-9093.266739736604*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-866.0254037844386*rdxCp2[0]*phiUx[1])-4330.127018922193*rdxCp2[0]*phiLy[1]+(750.0*phiUx[0]+750.0*phiLy[0]-3000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-6235.382907247957*rdxCp2Sq[0]*phiUx[1])-21823.84017536785*rdxCp2Sq[0]*phiLy[1]+(5400.0*phiUx[0]-10800.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-6235.382907247957*rdxCp2R3[0]*phiUx[1])-8183.940065762942*rdxCp2R3[0]*phiLy[1]+(5400.0*phiUx[0]-7875.0*phiLy[0]+4950.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[3])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2182.384017536785*rdxCp2R3[1])-6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]-2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-26292.53125889555*rdxCp2[0]*rdxCp2Sq[1])-12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]-415.6921938165305*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(75600.0*rdxCp2R4[1]+219960.0*rdxCp2[0]*rdxCp2R3[1]+89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-2182.384017536785*rdxCp2[0]*rdxCp2R3[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(32735.76026305177*rdxCp2R4[1]+94154.28189944413*rdxCp2[0]*rdxCp2R3[1]+37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(28350.0*phiLy[0]-66150.0*phiC[0])*rdxCp2R4[1]+((-36373.06695894642*rdxCp2[0]*phiUx[1])-24317.99333826703*rdxCp2[0]*phiLy[1]+(31500.0*phiUx[0]+81540.0*phiLy[0]-223020.0*phiC[0]-21000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-61279.95757178687*rdxCp2Sq[0]*phiUx[1])-11223.68923304632*rdxCp2Sq[0]*phiLy[1]+(53070.0*phiUx[0]+32355.0*phiLy[0]-130335.0*phiC[0]-111280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-21408.14798155132*rdxCp2R3[0]*phiUx[1])-311.7691453623978*rdxCp2R3[0]*phiLy[1]+(18540.0*phiUx[0]+900.0*phiLy[0]-21240.0*phiC[0]-48360.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-692.8203230275509*rdxCp2R4[0]*phiUx[1]+(600.0*phiUx[0]-600.0*phiC[0]-1600.0*bcVals[0])*rdxCp2R4[0])*omega+66150.0*phiC[0]*rdxCp2R4[1]+223020.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+600.0*phiC[0]*rdxCp2R4[0])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1770.0*rdxCp2[0]*rdxCp2Sq[1])-720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-31500.0*rdxCp2R3[1])-30300.0*rdxCp2[0]*rdxCp2Sq[1]-7740.0*rdxCp2Sq[0]*rdxCp2[1]-240.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+((346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(50402.67850025433*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3150.0*rdxCp2[0]*rdxCp2R3[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(23400.0*rdxCp2[0]*rdxCp2R3[1]+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(66150.0*phiC[1]-28350.0*phiLy[1])*rdxCp2R4[1]+(52500.0*rdxCp2[0]*phiUx[1]-25380.0*rdxCp2[0]*phiLy[1]+223020.0*rdxCp2[0]*phiC[1]+((-45466.33369868303*phiUx[0])+20264.99444855586*phiLy[0]-60621.7782649107*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(25200.0*rdxCp2Sq[0]*phiUx[1]-6435.0*rdxCp2Sq[0]*phiLy[1]+130335.0*rdxCp2Sq[0]*phiC[1]+((-21823.84017536785*phiUx[0])+9353.074360871933*phiLy[0]-72919.33899864974*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(900.0*rdxCp2R3[0]*phiUx[1]-180.0*rdxCp2R3[0]*phiLy[1]+21240.0*rdxCp2R3[0]*phiC[1]+((-779.4228634059946*phiUx[0])+259.8076211353315*phiLy[0]-21823.84017536785*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+600.0*rdxCp2R4[0]*phiC[1]-692.8203230275509*bcVals[0]*rdxCp2R4[0])*omega-66150.0*phiC[1]*rdxCp2R4[1]-223020.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-600.0*rdxCp2R4[0]*phiC[1]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+5455.960043841962*rho[0]*rdxCp2R3[1]+15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiUx[3]+((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-32735.76026305177*rdxCp2R4[1])-113795.7380572752*rdxCp2[0]*rdxCp2R3[1]-64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-5670.0*rdxCp2[0]*rdxCp2R3[1])-16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiUx[2]+(18900.0*rdxCp2R4[1]+59400.0*rdxCp2[0]*rdxCp2R3[1]+28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[2]+16367.88013152588*phiLy[0]*rdxCp2R4[1]+((-6300.0*rdxCp2[0]*phiUx[1])-6480.0*rdxCp2[0]*phiLy[1]+(5455.960043841962*phiUx[0]+51441.90898479563*phiLy[0]-3637.306695894642*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-8850.0*rdxCp2Sq[0]*phiUx[1])-6480.0*rdxCp2Sq[0]*phiLy[1]+(7664.324823492281*phiUx[0]+24811.62781842416*phiLy[0]-20438.19952931275*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3000.0*rdxCp2R3[0]*phiUx[1])-900.0*rdxCp2R3[0]*phiLy[1]+(2598.076211353316*phiUx[0]+2598.076211353316*phiLy[0]-7967.433714816835*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[2]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[3]+(4950.0*rdxCp2[0]*rdxCp2R3[1]-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(8183.940065762942*rdxCp2[0]*rdxCp2R3[1]+21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]+4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(6235.382907247957*rdxCp2[0]*rdxCp2R3[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-16367.88013152588*phiLy[1]*rdxCp2R4[1]+(9093.266739736604*rdxCp2[0]*phiUx[1]-36476.99000740053*rdxCp2[0]*phiLy[1]+((-7875.0*phiUx[0])+5400.0*phiLy[0]-10500.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((5400.0*phiLy[0]-8850.0*bcVals[0])*rdxCp2Sq[0]-9846.708841029065*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-866.0254037844386*rdxCp2R3[0]*phiUx[1])-519.6152422706631*rdxCp2R3[0]*phiLy[1]+(750.0*phiUx[0]+750.0*phiLy[0]-2600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[3])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-51.96152422706631*rdxCp2Sq[1])-571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+(571.5767664977294*rdxCp2[0]*rdxCp2[1]+51.96152422706631*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-200.0*rdxCp2R3[1])-2220.0*rdxCp2[0]*rdxCp2Sq[1]-100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-51.96152422706631*rdxCp2[0]*rdxCp2Sq[1])-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-86.60254037844386*rdxCp2R3[1])-987.26896031426*rdxCp2[0]*rdxCp2Sq[1]-173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(75.0*phiC[0]-75.0*phiLy[0])*rdxCp2R3[1]+(173.2050807568877*rdxCp2[0]*phiUx[1]+259.8076211353315*rdxCp2[0]*phiLy[1]+((-150.0*phiUx[0])-855.0*phiLy[0]+1005.0*phiC[0]+100.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(987.26896031426*rdxCp2Sq[0]*phiUx[1]+51.96152422706631*rdxCp2Sq[0]*phiLy[1]+((-855.0*phiUx[0])-150.0*phiLy[0]+1005.0*phiC[0]+2220.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+86.60254037844386*rdxCp2R3[0]*phiUx[1]+((-75.0*phiUx[0])+75.0*phiC[0]+200.0*bcVals[0])*rdxCp2R3[0])*omega-75.0*phiC[0]*rdxCp2R3[1]-1005.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-1005.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-75.0*phiC[0]*rdxCp2R3[0]))/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = (((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]-129.9038105676658*rdxCp2Sq[0]*rho[0])*omega*volFac+((259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1991.858428704209*rdxCp2[0]*rdxCp2Sq[1])-86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(225.0*rdxCp2[0]*rdxCp2Sq[1]-225.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-750.0*rdxCp2[0]*rdxCp2Sq[1])-150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(225.0*phiLy[1]-225.0*phiC[1])*rdxCp2R3[1]+((-750.0*rdxCp2[0]*phiUx[1])+765.0*rdxCp2[0]*phiLy[1]-3015.0*rdxCp2[0]*phiC[1]+(649.5190528383289*phiUx[0]-649.5190528383289*phiLy[0]+866.0254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-150.0*rdxCp2Sq[0]*phiUx[1])+90.0*rdxCp2Sq[0]*phiLy[1]-3015.0*rdxCp2Sq[0]*phiC[1]+(129.9038105676658*phiUx[0]-129.9038105676658*phiLy[0]+3031.088913245535*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-225.0*rdxCp2R3[0]*phiC[1]+259.8076211353315*bcVals[0]*rdxCp2R3[0])*omega+225.0*phiC[1]*rdxCp2R3[1]+3015.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+225.0*rdxCp2R3[0]*phiC[1])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-129.9038105676658*rho[0]*rdxCp2Sq[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiUx[3]+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-259.8076211353315*rdxCp2R3[1])-3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-90.0*rdxCp2[0]*rdxCp2Sq[1])-765.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiUx[2]+(150.0*rdxCp2[0]*rdxCp2Sq[1]+750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiC[2]+(150.0*rdxCp2[0]*phiUx[1]+225.0*rdxCp2[0]*phiLy[1]+((-129.9038105676658*phiUx[0])+129.9038105676658*phiLy[0]+86.60254037844386*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(750.0*rdxCp2Sq[0]*phiUx[1]-225.0*rdxCp2Sq[0]*phiLy[1]+((-649.5190528383289*phiUx[0])+649.5190528383289*phiLy[0]+1991.858428704209*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-225.0*rdxCp2R3[1])-3015.0*rdxCp2[0]*rdxCp2Sq[1]-3015.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiC[2]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+((-363.7306695894642*rdxCp2[0]*rdxCp2[1])-259.8076211353315*rdxCp2Sq[0])*rho[2]+(259.8076211353315*rdxCp2Sq[1]+363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-450.0*rdxCp2R3[1])-6030.0*rdxCp2[0]*rdxCp2Sq[1]-6030.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2R3[0])*phiC[3]+((-1300.0*rdxCp2[0]*rdxCp2Sq[1])-500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]+1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(433.0127018922193*rdxCp2Sq[0]*rdxCp2[1]-433.0127018922193*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-433.0127018922193*rdxCp2[0]*phiUx[1])-1299.038105676658*rdxCp2[0]*phiLy[1]+(375.0*phiUx[0]-375.0*phiLy[0]+500.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(433.0127018922193*rdxCp2Sq[0]*phiUx[1]-259.8076211353315*rdxCp2Sq[0]*phiLy[1]+((-375.0*phiUx[0])+375.0*phiLy[0]+1300.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0])*phiC[3])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-1454.922678357857*rdxCp2Sq[1])-8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+(8764.17708629852*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(21823.84017536785*rdxCp2R3[1]+134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]+18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-1091.192008768392*rdxCp2[0]*rdxCp2Sq[1])-8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-50400.0*rdxCp2R3[1])-314940.0*rdxCp2[0]*rdxCp2Sq[1]-63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(44100.0*phiC[0]-18900.0*phiUy[0])*rdxCp2R3[1]+(8105.997779422343*rdxCp2[0]*phiUy[1]-18186.53347947321*rdxCp2[0]*phiLx[1]-63000.0*rdxCp2[0]*bcVals[1]+((-116685.0*phiUy[0])-15750.0*phiLx[0]+321405.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(1091.192008768392*rdxCp2Sq[0]*phiUy[1]-134736.2323207829*rdxCp2Sq[0]*phiLx[1]-314940.0*rdxCp2Sq[0]*bcVals[1]+((-15750.0*phiUy[0])-116685.0*phiLx[0]+321405.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-21823.84017536785*rdxCp2R3[0]*phiLx[1]-50400.0*rdxCp2R3[0]*bcVals[1]+(44100.0*phiC[0]-18900.0*phiLx[0])*rdxCp2R3[0])*omega-44100.0*phiC[0]*rdxCp2R3[1]-321405.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-321405.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-44100.0*phiC[0]*rdxCp2R3[0]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = (((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]-3637.306695894642*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(23400.0*rdxCp2[0]*rdxCp2Sq[1]+3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-50402.67850025433*rdxCp2[0]*rdxCp2Sq[1])-10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(18900.0*phiUy[1]-44100.0*phiC[1])*rdxCp2R3[1]+(20745.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiLx[1]-321405.0*rdxCp2[0]*phiC[1]+181865.3347947321*rdxCp2[0]*bcVals[1]+((-20264.99444855586*phiUy[0])-45466.33369868303*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(2835.0*rdxCp2Sq[0]*phiUy[1]-85350.0*rdxCp2Sq[0]*phiLx[1]-321405.0*rdxCp2Sq[0]*phiC[1]+164198.4165575295*rdxCp2Sq[0]*bcVals[1]+((-2727.980021920981*phiUy[0])-73915.26821300182*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-12600.0*rdxCp2R3[0]*phiLx[1]-44100.0*rdxCp2R3[0]*phiC[1]+21823.84017536785*rdxCp2R3[0]*bcVals[1]-10911.92008768392*phiLx[0]*rdxCp2R3[0])*omega+44100.0*phiC[1]*rdxCp2R3[1]+321405.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+44100.0*rdxCp2R3[0]*phiC[1])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-3637.306695894642*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiLx[3]+(12600.0*rdxCp2R3[1]+85350.0*rdxCp2[0]*rdxCp2Sq[1]+52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-2835.0*rdxCp2[0]*rdxCp2Sq[1])-20745.0*rdxCp2Sq[0]*rdxCp2[1]-18900.0*rdxCp2R3[0])*phiLx[2]+(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0])*phiC[2]+(21823.84017536785*rdxCp2R3[1]+164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]+181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-10911.92008768392*phiUy[0]*rdxCp2R3[1]+(2160.0*rdxCp2[0]*phiUy[1]-3150.0*rdxCp2[0]*phiLx[1]-10911.92008768392*rdxCp2[0]*bcVals[1]+((-73915.26821300182*phiUy[0])-2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(3150.0*rdxCp2Sq[0]*phiUy[1]-23400.0*rdxCp2Sq[0]*phiLx[1]-50402.67850025433*rdxCp2Sq[0]*bcVals[1]+((-45466.33369868303*phiUy[0])-20264.99444855586*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-44100.0*rdxCp2R3[1])-321405.0*rdxCp2[0]*rdxCp2Sq[1]-321405.0*rdxCp2Sq[0]*rdxCp2[1]-44100.0*rdxCp2R3[0])*phiC[2]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-320.4293994002423*rdxCp2[0]*rdxCp2[1])-1212.435565298214*rdxCp2Sq[0])*rho[2]+(1212.435565298214*rdxCp2Sq[1]+320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]-1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiLx[3]+((-14700.0*rdxCp2R3[1])-107135.0*rdxCp2[0]*rdxCp2Sq[1]-107135.0*rdxCp2Sq[0]*rdxCp2[1]-14700.0*rdxCp2R3[0])*phiC[3]+(2078.460969082652*rdxCp2[0]*rdxCp2Sq[1]+3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-2727.980021920981*rdxCp2[0]*rdxCp2Sq[1])-16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]-3637.306695894642*rdxCp2R3[0])*phiLx[2]+(10500.0*rdxCp2Sq[0]*rdxCp2[1]-1650.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+3637.306695894642*phiUy[1]*rdxCp2R3[1]+(16116.7327644284*rdxCp2[0]*phiUy[1]-3031.088913245535*rdxCp2[0]*phiLx[1]+10500.0*rdxCp2[0]*bcVals[1]+((-1800.0*phiUy[0])-2625.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(2727.980021920981*rdxCp2Sq[0]*phiUy[1]-2078.460969082652*rdxCp2Sq[0]*phiLx[1]-1650.0*rdxCp2Sq[0]*bcVals[1]+((-2625.0*phiUy[0])-1800.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0])*phiC[3])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-415.6921938165305*rdxCp2R3[1])-12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]-26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-2286.307065990918*rdxCp2[0]*rdxCp2Sq[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]-2182.384017536785*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+((1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-692.8203230275509*rdxCp2R4[1])-21408.14798155132*rdxCp2[0]*rdxCp2R3[1]-61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2R3[1])-11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]-24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-1600.0*rdxCp2R4[1])-48360.0*rdxCp2[0]*rdxCp2R3[1]-111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(600.0*phiUy[0]-600.0*phiC[0])*rdxCp2R4[1]+((-1039.230484541326*rdxCp2[0]*phiUy[1])+1039.230484541326*rdxCp2[0]*phiLx[1]+3600.0*rdxCp2[0]*bcVals[1]+(18540.0*phiUy[0]+900.0*phiLx[0]-21240.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-3065.729929396912*rdxCp2Sq[0]*phiUy[1])+37360.33591926068*rdxCp2Sq[0]*phiLx[1]+89820.0*rdxCp2Sq[0]*bcVals[1]+(53070.0*phiUy[0]+32355.0*phiLx[0]-130335.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2182.384017536785*rdxCp2R3[0]*phiUy[1])+94154.28189944413*rdxCp2R3[0]*phiLx[1]+219960.0*rdxCp2R3[0]*bcVals[1]+(31500.0*phiUy[0]+81540.0*phiLx[0]-223020.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+32735.76026305177*rdxCp2R4[0]*phiLx[1]+75600.0*rdxCp2R4[0]*bcVals[1]+(28350.0*phiLx[0]-66150.0*phiC[0])*rdxCp2R4[0])*omega+600.0*phiC[0]*rdxCp2R4[1]+21240.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+66150.0*phiC[0]*rdxCp2R4[0])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]+5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0]*rho[0])*omega*volFac+((692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3000.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6480.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-7967.433714816835*rdxCp2[0]*rdxCp2R3[1])-20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(600.0*phiC[1]-600.0*phiUy[1])*rdxCp2R4[1]+((-6240.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiLx[1]+21240.0*rdxCp2[0]*phiC[1]-10392.30484541326*rdxCp2[0]*bcVals[1]+(2598.076211353316*phiUy[0]+2598.076211353316*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-16785.0*rdxCp2Sq[0]*phiUy[1])+28650.0*rdxCp2Sq[0]*phiLx[1]+130335.0*rdxCp2Sq[0]*phiC[1]-64951.90528383289*rdxCp2Sq[0]*bcVals[1]+(7664.324823492281*phiUy[0]+24811.62781842416*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-5670.0*rdxCp2R3[0]*phiUy[1])+59400.0*rdxCp2R3[0]*phiLx[1]+223020.0*rdxCp2R3[0]*phiC[1]-113795.7380572752*rdxCp2R3[0]*bcVals[1]+(5455.960043841962*phiUy[0]+51441.90898479563*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+18900.0*rdxCp2R4[0]*phiLx[1]+66150.0*rdxCp2R4[0]*phiC[1]-32735.76026305177*rdxCp2R4[0]*bcVals[1]+16367.88013152588*phiLx[0]*rdxCp2R4[0])*omega-600.0*phiC[1]*rdxCp2R4[1]-21240.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-66150.0*rdxCp2R4[0]*phiC[1]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+((-240.0*rdxCp2R3[1])-7740.0*rdxCp2[0]*rdxCp2Sq[1]-30300.0*rdxCp2Sq[0]*rdxCp2[1]-31500.0*rdxCp2R3[0])*rho[2]+((-720.0*rdxCp2[0]*rdxCp2Sq[1])-1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiLx[3]+(900.0*rdxCp2[0]*rdxCp2R3[1]+25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+52500.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-180.0*rdxCp2[0]*rdxCp2R3[1])-6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25380.0*rdxCp2R3[0]*rdxCp2[1]-28350.0*rdxCp2R4[0])*phiLx[2]+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[2]+((-692.8203230275509*rdxCp2R4[1])-21823.84017536785*rdxCp2[0]*rdxCp2R3[1]-72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]-60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-300.0*rdxCp2[0]*phiUy[1])+300.0*rdxCp2[0]*phiLx[1]+1039.230484541326*rdxCp2[0]*bcVals[1]+(259.8076211353315*phiLx[0]-779.4228634059946*phiUy[0])*rdxCp2[0])*rdxCp2R3[1]+(10800.0*rdxCp2Sq[0]*phiLx[1]+24941.53162899183*rdxCp2Sq[0]*bcVals[1]+(9353.074360871933*phiLx[0]-21823.84017536785*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3150.0*rdxCp2R3[0]*phiUy[1]+23400.0*rdxCp2R3[0]*phiLx[1]+50402.67850025433*rdxCp2R3[0]*bcVals[1]+(20264.99444855586*phiLx[0]-45466.33369868303*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[2]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiLx[3]+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[3]+(9093.266739736604*rdxCp2R3[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2R3[1])-9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]-36476.99000740053*rdxCp2R3[0]*rdxCp2[1]-16367.88013152588*rdxCp2R4[0])*phiLx[2]+((-2600.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(4330.127018922193*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiLx[1]-3000.0*rdxCp2[0]*bcVals[1]+(750.0*phiUy[0]+750.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]+6235.382907247957*rdxCp2Sq[0]*phiLx[1]-10800.0*rdxCp2Sq[0]*bcVals[1]+5400.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(8183.940065762942*rdxCp2R3[0]*phiUy[1]+6235.382907247957*rdxCp2R3[0]*phiLx[1]+4950.0*rdxCp2R3[0]*bcVals[1]+(5400.0*phiLx[0]-7875.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[3])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = (((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2182.384017536785*rdxCp2R3[1]+6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]+2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(26292.53125889555*rdxCp2[0]*rdxCp2Sq[1]+12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-32735.76026305177*rdxCp2R4[1])-94154.28189944413*rdxCp2[0]*rdxCp2R3[1]-37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(2182.384017536785*rdxCp2[0]*rdxCp2R3[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(75600.0*rdxCp2R4[1]+219960.0*rdxCp2[0]*rdxCp2R3[1]+89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(28350.0*phiUy[0]-66150.0*phiC[0])*rdxCp2R4[1]+(24317.99333826703*rdxCp2[0]*phiUy[1]+36373.06695894642*rdxCp2[0]*phiLx[1]+21000.0*rdxCp2[0]*bcVals[1]+(81540.0*phiUy[0]+31500.0*phiLx[0]-223020.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(11223.68923304632*rdxCp2Sq[0]*phiUy[1]+61279.95757178687*rdxCp2Sq[0]*phiLx[1]+111280.0*rdxCp2Sq[0]*bcVals[1]+(32355.0*phiUy[0]+53070.0*phiLx[0]-130335.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(311.7691453623978*rdxCp2R3[0]*phiUy[1]+21408.14798155132*rdxCp2R3[0]*phiLx[1]+48360.0*rdxCp2R3[0]*bcVals[1]+(900.0*phiUy[0]+18540.0*phiLx[0]-21240.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+692.8203230275509*rdxCp2R4[0]*phiLx[1]+1600.0*rdxCp2R4[0]*bcVals[1]+(600.0*phiLx[0]-600.0*phiC[0])*rdxCp2R4[0])*omega+66150.0*phiC[0]*rdxCp2R4[1]+223020.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+600.0*phiC[0]*rdxCp2R4[0])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = (((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1770.0*rdxCp2[0]*rdxCp2Sq[1]+720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(31500.0*rdxCp2R3[1]+30300.0*rdxCp2[0]*rdxCp2Sq[1]+7740.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+(((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-23400.0*rdxCp2[0]*rdxCp2R3[1])-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(300.0*rdxCp2R3[0]*rdxCp2[1]-3150.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(50402.67850025433*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(28350.0*phiUy[1]-66150.0*phiC[1])*rdxCp2R4[1]+(25380.0*rdxCp2[0]*phiUy[1]-52500.0*rdxCp2[0]*phiLx[1]-223020.0*rdxCp2[0]*phiC[1]+60621.7782649107*rdxCp2[0]*bcVals[1]+(20264.99444855586*phiUy[0]-45466.33369868303*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(6435.0*rdxCp2Sq[0]*phiUy[1]-25200.0*rdxCp2Sq[0]*phiLx[1]-130335.0*rdxCp2Sq[0]*phiC[1]+72919.33899864974*rdxCp2Sq[0]*bcVals[1]+(9353.074360871933*phiUy[0]-21823.84017536785*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(180.0*rdxCp2R3[0]*phiUy[1]-900.0*rdxCp2R3[0]*phiLx[1]-21240.0*rdxCp2R3[0]*phiC[1]+21823.84017536785*rdxCp2R3[0]*bcVals[1]+(259.8076211353315*phiUy[0]-779.4228634059946*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-600.0*rdxCp2R4[0]*phiC[1]+692.8203230275509*rdxCp2R4[0]*bcVals[1])*omega+66150.0*phiC[1]*rdxCp2R4[1]+223020.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+600.0*rdxCp2R4[0]*phiC[1])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = (((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+5455.960043841962*rho[0]*rdxCp2R3[1]+15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiLx[3]+((-18900.0*rdxCp2R4[1])-59400.0*rdxCp2[0]*rdxCp2R3[1]-28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(5670.0*rdxCp2[0]*rdxCp2R3[1]+16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiLx[2]+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[2]+((-32735.76026305177*rdxCp2R4[1])-113795.7380572752*rdxCp2[0]*rdxCp2R3[1]-64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+16367.88013152588*phiUy[0]*rdxCp2R4[1]+(6480.0*rdxCp2[0]*phiUy[1]+6300.0*rdxCp2[0]*phiLx[1]+3637.306695894642*rdxCp2[0]*bcVals[1]+(51441.90898479563*phiUy[0]+5455.960043841962*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(6480.0*rdxCp2Sq[0]*phiUy[1]+8850.0*rdxCp2Sq[0]*phiLx[1]+20438.19952931275*rdxCp2Sq[0]*bcVals[1]+(24811.62781842416*phiUy[0]+7664.324823492281*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(900.0*rdxCp2R3[0]*phiUy[1]+3000.0*rdxCp2R3[0]*phiLx[1]+7967.433714816835*rdxCp2R3[0]*bcVals[1]+(2598.076211353316*phiUy[0]+2598.076211353316*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[2])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[3]+((-6235.382907247957*rdxCp2[0]*rdxCp2R3[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-8183.940065762942*rdxCp2[0]*rdxCp2R3[1])-21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]-4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(4950.0*rdxCp2[0]*rdxCp2R3[1]-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+16367.88013152588*phiUy[1]*rdxCp2R4[1]+(36476.99000740053*rdxCp2[0]*phiUy[1]-9093.266739736604*rdxCp2[0]*phiLx[1]+10500.0*rdxCp2[0]*bcVals[1]+(5400.0*phiUy[0]-7875.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(9846.708841029065*rdxCp2Sq[0]*phiUy[1]+8850.0*rdxCp2Sq[0]*bcVals[1]+5400.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2R3[0]*phiUy[1]+866.0254037844386*rdxCp2R3[0]*phiLx[1]+2600.0*rdxCp2R3[0]*bcVals[1]+(750.0*phiUy[0]+750.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[3])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(51.96152422706631*rdxCp2Sq[1]+571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+((-571.5767664977294*rdxCp2[0]*rdxCp2[1])-51.96152422706631*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(86.60254037844386*rdxCp2R3[1]+987.26896031426*rdxCp2[0]*rdxCp2Sq[1]+173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(51.96152422706631*rdxCp2[0]*rdxCp2Sq[1]+259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(200.0*rdxCp2R3[1]+2220.0*rdxCp2[0]*rdxCp2Sq[1]+100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(75.0*phiC[0]-75.0*phiUy[0])*rdxCp2R3[1]+((-259.8076211353315*rdxCp2[0]*phiUy[1])-173.2050807568877*rdxCp2[0]*phiLx[1]-100.0*rdxCp2[0]*bcVals[1]+((-855.0*phiUy[0])-150.0*phiLx[0]+1005.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-51.96152422706631*rdxCp2Sq[0]*phiUy[1])-987.26896031426*rdxCp2Sq[0]*phiLx[1]-2220.0*rdxCp2Sq[0]*bcVals[1]+((-150.0*phiUy[0])-855.0*phiLx[0]+1005.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-86.60254037844386*rdxCp2R3[0]*phiLx[1]-200.0*rdxCp2R3[0]*bcVals[1]+(75.0*phiC[0]-75.0*phiLx[0])*rdxCp2R3[0])*omega-75.0*phiC[0]*rdxCp2R3[1]-1005.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-1005.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-75.0*phiC[0]*rdxCp2R3[0]))/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]-129.9038105676658*rdxCp2Sq[0]*rho[0])*omega*volFac+((259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(750.0*rdxCp2[0]*rdxCp2Sq[1]+150.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(225.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(1991.858428704209*rdxCp2[0]*rdxCp2Sq[1]+86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(225.0*phiC[1]-225.0*phiUy[1])*rdxCp2R3[1]+((-765.0*rdxCp2[0]*phiUy[1])+750.0*rdxCp2[0]*phiLx[1]+3015.0*rdxCp2[0]*phiC[1]-866.0254037844386*rdxCp2[0]*bcVals[1]+(649.5190528383289*phiLx[0]-649.5190528383289*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-90.0*rdxCp2Sq[0]*phiUy[1])+150.0*rdxCp2Sq[0]*phiLx[1]+3015.0*rdxCp2Sq[0]*phiC[1]-3031.088913245535*rdxCp2Sq[0]*bcVals[1]+(129.9038105676658*phiLx[0]-129.9038105676658*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]+225.0*rdxCp2R3[0]*phiC[1]-259.8076211353315*rdxCp2R3[0]*bcVals[1])*omega-225.0*phiC[1]*rdxCp2R3[1]-3015.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-3015.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-225.0*rdxCp2R3[0]*phiC[1]))/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = (((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-129.9038105676658*rho[0]*rdxCp2Sq[1]-1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiLx[3]+((-150.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(90.0*rdxCp2[0]*rdxCp2Sq[1]+765.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiLx[2]+((-225.0*rdxCp2R3[1])-3015.0*rdxCp2[0]*rdxCp2Sq[1]-3015.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiC[2]+(259.8076211353315*rdxCp2R3[1]+3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-225.0*rdxCp2[0]*phiUy[1])-150.0*rdxCp2[0]*phiLx[1]-86.60254037844386*rdxCp2[0]*bcVals[1]+(129.9038105676658*phiUy[0]-129.9038105676658*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(225.0*rdxCp2Sq[0]*phiUy[1]-750.0*rdxCp2Sq[0]*phiLx[1]-1991.858428704209*rdxCp2Sq[0]*bcVals[1]+(649.5190528383289*phiUy[0]-649.5190528383289*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiC[2])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+(363.7306695894642*rdxCp2[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0])*rho[2]+((-259.8076211353315*rdxCp2Sq[1])-363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-450.0*rdxCp2R3[1])-6030.0*rdxCp2[0]*rdxCp2Sq[1]-6030.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2R3[0])*phiC[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])-1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(1300.0*rdxCp2[0]*rdxCp2Sq[1]+500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(1299.038105676658*rdxCp2[0]*phiUy[1]+433.0127018922193*rdxCp2[0]*phiLx[1]-500.0*rdxCp2[0]*bcVals[1]+(375.0*phiLx[0]-375.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(259.8076211353315*rdxCp2Sq[0]*phiUy[1]-433.0127018922193*rdxCp2Sq[0]*phiLx[1]-1300.0*rdxCp2Sq[0]*bcVals[1]+(375.0*phiUy[0]-375.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0])*phiC[3])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((708.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-1454.922678357857*rdxCp2Sq[1])-8764.17708629852*rdxCp2[0]*rdxCp2[1])*rho[2]+((-8764.17708629852*rdxCp2[0]*rdxCp2[1])-1454.922678357857*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-9360.0*rdxCp2[0]*rdxCp2Sq[1])-1260.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1260.0*rdxCp2[0]*rdxCp2Sq[1])-9360.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(50400.0*rdxCp2R3[1]+314940.0*rdxCp2[0]*rdxCp2Sq[1]+63000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(21823.84017536785*rdxCp2R3[1]+134736.2323207829*rdxCp2[0]*rdxCp2Sq[1]+18186.53347947321*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-1091.192008768392*rdxCp2[0]*rdxCp2Sq[1])-8105.997779422343*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(18900.0*phiLy[0]-44100.0*phiC[0])*rdxCp2R3[1]+((-8105.997779422343*rdxCp2[0]*phiLy[1])+18186.53347947321*rdxCp2[0]*phiLx[1]+63000.0*rdxCp2[0]*bcVals[1]+(116685.0*phiLy[0]+15750.0*phiLx[0]-321405.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-1091.192008768392*rdxCp2Sq[0]*phiLy[1])+134736.2323207829*rdxCp2Sq[0]*phiLx[1]+314940.0*rdxCp2Sq[0]*bcVals[1]+(15750.0*phiLy[0]+116685.0*phiLx[0]-321405.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+21823.84017536785*rdxCp2R3[0]*phiLx[1]+50400.0*rdxCp2R3[0]*bcVals[1]+(18900.0*phiLx[0]-44100.0*phiC[0])*rdxCp2R3[0])*omega+44100.0*phiC[0]*rdxCp2R3[1]+321405.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+321405.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+44100.0*phiC[0]*rdxCp2R3[0])/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((1454.922678357857*rdxCp2Sq[1]+384.5152792802907*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1]+3637.306695894642*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-21823.84017536785*rdxCp2R3[1])-23954.26266867757*rdxCp2[0]*rdxCp2Sq[1]-3273.576026305177*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-3637.306695894642*rdxCp2[0]*rdxCp2Sq[1])-2494.153162899182*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(50402.67850025433*rdxCp2[0]*rdxCp2Sq[1]+10911.92008768392*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(23400.0*rdxCp2[0]*rdxCp2Sq[1]+3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(44100.0*phiC[1]-18900.0*phiLy[1])*rdxCp2R3[1]+((-20745.0*rdxCp2[0]*phiLy[1])+52500.0*rdxCp2[0]*phiLx[1]+321405.0*rdxCp2[0]*phiC[1]-181865.3347947321*rdxCp2[0]*bcVals[1]+(20264.99444855586*phiLy[0]+45466.33369868303*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-2835.0*rdxCp2Sq[0]*phiLy[1])+85350.0*rdxCp2Sq[0]*phiLx[1]+321405.0*rdxCp2Sq[0]*phiC[1]-164198.4165575295*rdxCp2Sq[0]*bcVals[1]+(2727.980021920981*phiLy[0]+73915.26821300182*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+12600.0*rdxCp2R3[0]*phiLx[1]+44100.0*rdxCp2R3[0]*phiC[1]-21823.84017536785*rdxCp2R3[0]*bcVals[1]+10911.92008768392*phiLx[0]*rdxCp2R3[0])*omega-44100.0*phiC[1]*rdxCp2R3[1]-321405.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-321405.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-44100.0*rdxCp2R3[0]*phiC[1]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((384.5152792802907*rdxCp2[0]*rdxCp2[1]+1454.922678357857*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+3637.306695894642*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-2494.153162899182*rdxCp2[0]*rdxCp2Sq[1])-3637.306695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-3273.576026305177*rdxCp2[0]*rdxCp2Sq[1])-23954.26266867757*rdxCp2Sq[0]*rdxCp2[1]-21823.84017536785*rdxCp2R3[0])*phiLx[3]+((-21823.84017536785*rdxCp2R3[1])-164198.4165575295*rdxCp2[0]*rdxCp2Sq[1]-181865.3347947321*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(12600.0*rdxCp2R3[1]+85350.0*rdxCp2[0]*rdxCp2Sq[1]+52500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-2835.0*rdxCp2[0]*rdxCp2Sq[1])-20745.0*rdxCp2Sq[0]*rdxCp2[1]-18900.0*rdxCp2R3[0])*phiLx[2]+(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0])*phiC[2]+10911.92008768392*phiLy[0]*rdxCp2R3[1]+((-2160.0*rdxCp2[0]*phiLy[1])+3150.0*rdxCp2[0]*phiLx[1]+10911.92008768392*rdxCp2[0]*bcVals[1]+(73915.26821300182*phiLy[0]+2727.980021920981*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-3150.0*rdxCp2Sq[0]*phiLy[1])+23400.0*rdxCp2Sq[0]*phiLx[1]+50402.67850025433*rdxCp2Sq[0]*bcVals[1]+(45466.33369868303*phiLy[0]+20264.99444855586*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-44100.0*rdxCp2R3[1])-321405.0*rdxCp2[0]*rdxCp2Sq[1]-321405.0*rdxCp2Sq[0]*rdxCp2[1]-44100.0*rdxCp2R3[0])*phiC[2]))/(44100.0*rdxCp2R3[1]+321405.0*rdxCp2[0]*rdxCp2Sq[1]+321405.0*rdxCp2Sq[0]*rdxCp2[1]+44100.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-320.4293994002423*rdxCp2[0]*rdxCp2[1])-1212.435565298214*rdxCp2Sq[0])*rho[2]+((-1212.435565298214*rdxCp2Sq[1])-320.4293994002423*rdxCp2[0]*rdxCp2[1])*rho[1]+1475.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-4200.0*rdxCp2R3[1])-18610.0*rdxCp2[0]*rdxCp2Sq[1]-3150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-3150.0*rdxCp2[0]*rdxCp2Sq[1])-18610.0*rdxCp2Sq[0]*rdxCp2[1]-4200.0*rdxCp2R3[0])*phiLx[3]+((-14700.0*rdxCp2R3[1])-107135.0*rdxCp2[0]*rdxCp2Sq[1]-107135.0*rdxCp2Sq[0]*rdxCp2[1]-14700.0*rdxCp2R3[0])*phiC[3]+(1650.0*rdxCp2[0]*rdxCp2Sq[1]-10500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2078.460969082652*rdxCp2[0]*rdxCp2Sq[1]+3031.088913245535*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-2727.980021920981*rdxCp2[0]*rdxCp2Sq[1])-16116.7327644284*rdxCp2Sq[0]*rdxCp2[1]-3637.306695894642*rdxCp2R3[0])*phiLx[2]-3637.306695894642*phiLy[1]*rdxCp2R3[1]+((-16116.7327644284*rdxCp2[0]*phiLy[1])+3031.088913245535*rdxCp2[0]*phiLx[1]-10500.0*rdxCp2[0]*bcVals[1]+(1800.0*phiLy[0]+2625.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-2727.980021920981*rdxCp2Sq[0]*phiLy[1])+2078.460969082652*rdxCp2Sq[0]*phiLx[1]+1650.0*rdxCp2Sq[0]*bcVals[1]+(2625.0*phiLy[0]+1800.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0])*phiC[3])/(14700.0*rdxCp2R3[1]+107135.0*rdxCp2[0]*rdxCp2Sq[1]+107135.0*rdxCp2Sq[0]*rdxCp2[1]+14700.0*rdxCp2R3[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((864.0*rdxCp2[0]*rdxCp2Sq[1]+2124.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-415.6921938165305*rdxCp2R3[1])-12470.76581449591*rdxCp2[0]*rdxCp2Sq[1]-26292.53125889555*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(2286.307065990918*rdxCp2[0]*rdxCp2Sq[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+((1200.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2520.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-28080.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-1600.0*rdxCp2R4[1])-48360.0*rdxCp2[0]*rdxCp2R3[1]-111280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-692.8203230275509*rdxCp2R4[1])-21408.14798155132*rdxCp2[0]*rdxCp2R3[1]-61279.95757178687*rdxCp2Sq[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-311.7691453623978*rdxCp2[0]*rdxCp2R3[1])-11223.68923304632*rdxCp2Sq[0]*rdxCp2Sq[1]-24317.99333826703*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(600.0*phiC[0]-600.0*phiLy[0])*rdxCp2R4[1]+(1039.230484541326*rdxCp2[0]*phiLy[1]-1039.230484541326*rdxCp2[0]*phiLx[1]-3600.0*rdxCp2[0]*bcVals[1]+((-18540.0*phiLy[0])-900.0*phiLx[0]+21240.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(3065.729929396912*rdxCp2Sq[0]*phiLy[1]-37360.33591926068*rdxCp2Sq[0]*phiLx[1]-89820.0*rdxCp2Sq[0]*bcVals[1]+((-53070.0*phiLy[0])-32355.0*phiLx[0]+130335.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(2182.384017536785*rdxCp2R3[0]*phiLy[1]-94154.28189944413*rdxCp2R3[0]*phiLx[1]-219960.0*rdxCp2R3[0]*bcVals[1]+((-31500.0*phiLy[0])-81540.0*phiLx[0]+223020.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-32735.76026305177*rdxCp2R4[0]*phiLx[1]-75600.0*rdxCp2R4[0]*bcVals[1]+(66150.0*phiC[0]-28350.0*phiLx[0])*rdxCp2R4[0])*omega-600.0*phiC[0]*rdxCp2R4[1]-21240.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-66150.0*phiC[0]*rdxCp2R4[0]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[1] = (((415.6921938165305*rdxCp2R3[1]+2244.737846609264*rdxCp2[0]*rdxCp2Sq[1]+1153.545837840872*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]-5715.767664977294*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-15328.64964698456*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0]*rho[0])*omega*volFac+((692.8203230275509*rdxCp2R4[1]+7205.331359486529*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+6547.152052610353*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-1039.230484541326*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-7482.459488697546*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-7967.433714816835*rdxCp2[0]*rdxCp2R3[1])-20438.19952931275*rdxCp2Sq[0]*rdxCp2Sq[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-3000.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-900.0*rdxCp2[0]*rdxCp2R3[1])-6480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-6480.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(600.0*phiLy[1]-600.0*phiC[1])*rdxCp2R4[1]+(6240.0*rdxCp2[0]*phiLy[1]-3000.0*rdxCp2[0]*phiLx[1]-21240.0*rdxCp2[0]*phiC[1]+10392.30484541326*rdxCp2[0]*bcVals[1]+((-2598.076211353316*phiLy[0])-2598.076211353316*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(16785.0*rdxCp2Sq[0]*phiLy[1]-28650.0*rdxCp2Sq[0]*phiLx[1]-130335.0*rdxCp2Sq[0]*phiC[1]+64951.90528383289*rdxCp2Sq[0]*bcVals[1]+((-7664.324823492281*phiLy[0])-24811.62781842416*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(5670.0*rdxCp2R3[0]*phiLy[1]-59400.0*rdxCp2R3[0]*phiLx[1]-223020.0*rdxCp2R3[0]*phiC[1]+113795.7380572752*rdxCp2R3[0]*bcVals[1]+((-5455.960043841962*phiLy[0])-51441.90898479563*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-18900.0*rdxCp2R4[0]*phiLx[1]-66150.0*rdxCp2R4[0]*phiC[1]+32735.76026305177*rdxCp2R4[0]*bcVals[1]-16367.88013152588*phiLx[0]*rdxCp2R4[0])*omega+600.0*phiC[1]*rdxCp2R4[1]+21240.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+66150.0*rdxCp2R4[0]*phiC[1])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((290.9845356715713*rdxCp2[0]*rdxCp2Sq[1]+1226.291971758765*rdxCp2Sq[0]*rdxCp2[1]+2182.384017536785*rdxCp2R3[0])*rho[3]+((-240.0*rdxCp2R3[1])-7740.0*rdxCp2[0]*rdxCp2Sq[1]-30300.0*rdxCp2Sq[0]*rdxCp2[1]-31500.0*rdxCp2R3[0])*rho[2]+(720.0*rdxCp2[0]*rdxCp2Sq[1]+1770.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((346.4101615137754*rdxCp2[0]*rdxCp2R3[1]-3637.306695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-207.8460969082653*rdxCp2[0]*rdxCp2R3[1])-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-29306.29966406539*rdxCp2R3[0]*rdxCp2[1]-32735.76026305177*rdxCp2R4[0])*phiLx[3]+((-692.8203230275509*rdxCp2R4[1])-21823.84017536785*rdxCp2[0]*rdxCp2R3[1]-72919.33899864974*rdxCp2Sq[0]*rdxCp2Sq[1]-60621.7782649107*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(900.0*rdxCp2[0]*rdxCp2R3[1]+25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]+52500.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-180.0*rdxCp2[0]*rdxCp2R3[1])-6435.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25380.0*rdxCp2R3[0]*rdxCp2[1]-28350.0*rdxCp2R4[0])*phiLx[2]+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[2]+(300.0*rdxCp2[0]*phiLy[1]-300.0*rdxCp2[0]*phiLx[1]-1039.230484541326*rdxCp2[0]*bcVals[1]+(779.4228634059946*phiLy[0]-259.8076211353315*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-10800.0*rdxCp2Sq[0]*phiLx[1])-24941.53162899183*rdxCp2Sq[0]*bcVals[1]+(21823.84017536785*phiLy[0]-9353.074360871933*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3150.0*rdxCp2R3[0]*phiLy[1])-23400.0*rdxCp2R3[0]*phiLx[1]-50402.67850025433*rdxCp2R3[0]*bcVals[1]+(45466.33369868303*phiLy[0]-20264.99444855586*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[2]))/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 
  phiC[3] = (((240.0*rdxCp2R3[1]+4296.0*rdxCp2[0]*rdxCp2Sq[1]+15786.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-1800.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-4425.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-5000.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-9450.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-600.0*rdxCp2[0]*rdxCp2R3[1])-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42120.0*rdxCp2R3[0]*rdxCp2[1]-18900.0*rdxCp2R4[0])*phiLx[3]+((-600.0*rdxCp2R4[1])-21240.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-223020.0*rdxCp2R3[0]*rdxCp2[1]-66150.0*rdxCp2R4[0])*phiC[3]+((-2600.0*rdxCp2[0]*rdxCp2R3[1])-8850.0*rdxCp2Sq[0]*rdxCp2Sq[1]-10500.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(9093.266739736604*rdxCp2R3[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2R3[1])-9846.708841029065*rdxCp2Sq[0]*rdxCp2Sq[1]-36476.99000740053*rdxCp2R3[0]*rdxCp2[1]-16367.88013152588*rdxCp2R4[0])*phiLx[2]+((-4330.127018922193*rdxCp2[0]*phiLy[1])-866.0254037844386*rdxCp2[0]*phiLx[1]+3000.0*rdxCp2[0]*bcVals[1]+((-750.0*phiLy[0])-750.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-21823.84017536785*rdxCp2Sq[0]*phiLy[1])-6235.382907247957*rdxCp2Sq[0]*phiLx[1]+10800.0*rdxCp2Sq[0]*bcVals[1]-5400.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-8183.940065762942*rdxCp2R3[0]*phiLy[1])-6235.382907247957*rdxCp2R3[0]*phiLx[1]-4950.0*rdxCp2R3[0]*bcVals[1]+(7875.0*phiLy[0]-5400.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0])*phiC[3])/(600.0*rdxCp2R4[1]+21240.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+223020.0*rdxCp2R3[0]*rdxCp2[1]+66150.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = -(1.0*(((2124.0*rdxCp2[0]*rdxCp2Sq[1]+864.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2182.384017536785*rdxCp2R3[1]+6131.459858793824*rdxCp2[0]*rdxCp2Sq[1]+2286.307065990918*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-26292.53125889555*rdxCp2[0]*rdxCp2Sq[1])-12470.76581449591*rdxCp2Sq[0]*rdxCp2[1]-415.6921938165305*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-28080.0*rdxCp2[0]*rdxCp2R3[1])-12960.0*rdxCp2Sq[0]*rdxCp2Sq[1]-360.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(2520.0*rdxCp2[0]*rdxCp2R3[1]+3540.0*rdxCp2Sq[0]*rdxCp2Sq[1]+1200.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-75600.0*rdxCp2R4[1])-219960.0*rdxCp2[0]*rdxCp2R3[1]-89820.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3600.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-32735.76026305177*rdxCp2R4[1])-94154.28189944413*rdxCp2[0]*rdxCp2R3[1]-37360.33591926068*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(2182.384017536785*rdxCp2[0]*rdxCp2R3[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2Sq[1]+1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(66150.0*phiC[0]-28350.0*phiLy[0])*rdxCp2R4[1]+((-24317.99333826703*rdxCp2[0]*phiLy[1])-36373.06695894642*rdxCp2[0]*phiLx[1]-21000.0*rdxCp2[0]*bcVals[1]+((-81540.0*phiLy[0])-31500.0*phiLx[0]+223020.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-11223.68923304632*rdxCp2Sq[0]*phiLy[1])-61279.95757178687*rdxCp2Sq[0]*phiLx[1]-111280.0*rdxCp2Sq[0]*bcVals[1]+((-32355.0*phiLy[0])-53070.0*phiLx[0]+130335.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-311.7691453623978*rdxCp2R3[0]*phiLy[1])-21408.14798155132*rdxCp2R3[0]*phiLx[1]-48360.0*rdxCp2R3[0]*bcVals[1]+((-900.0*phiLy[0])-18540.0*phiLx[0]+21240.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-692.8203230275509*rdxCp2R4[0]*phiLx[1]-1600.0*rdxCp2R4[0]*bcVals[1]+(600.0*phiC[0]-600.0*phiLx[0])*rdxCp2R4[0])*omega-66150.0*phiC[0]*rdxCp2R4[1]-223020.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-130335.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-600.0*phiC[0]*rdxCp2R4[0]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((2182.384017536785*rdxCp2R3[1]+1226.291971758765*rdxCp2[0]*rdxCp2Sq[1]+290.9845356715713*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1770.0*rdxCp2[0]*rdxCp2Sq[1]+720.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-31500.0*rdxCp2R3[1])-30300.0*rdxCp2[0]*rdxCp2Sq[1]-7740.0*rdxCp2Sq[0]*rdxCp2[1]-240.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+(((-32735.76026305177*rdxCp2R4[1])-29306.29966406539*rdxCp2[0]*rdxCp2R3[1]-7430.497964470483*rdxCp2Sq[0]*rdxCp2Sq[1]-207.8460969082653*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(346.4101615137754*rdxCp2R3[0]*rdxCp2[1]-3637.306695894642*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-50402.67850025433*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-23400.0*rdxCp2[0]*rdxCp2R3[1])-10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-300.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(300.0*rdxCp2R3[0]*rdxCp2[1]-3150.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(66150.0*phiC[1]-28350.0*phiLy[1])*rdxCp2R4[1]+((-25380.0*rdxCp2[0]*phiLy[1])+52500.0*rdxCp2[0]*phiLx[1]+223020.0*rdxCp2[0]*phiC[1]-60621.7782649107*rdxCp2[0]*bcVals[1]+(45466.33369868303*phiLx[0]-20264.99444855586*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-6435.0*rdxCp2Sq[0]*phiLy[1])+25200.0*rdxCp2Sq[0]*phiLx[1]+130335.0*rdxCp2Sq[0]*phiC[1]-72919.33899864974*rdxCp2Sq[0]*bcVals[1]+(21823.84017536785*phiLx[0]-9353.074360871933*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-180.0*rdxCp2R3[0]*phiLy[1])+900.0*rdxCp2R3[0]*phiLx[1]+21240.0*rdxCp2R3[0]*phiC[1]-21823.84017536785*rdxCp2R3[0]*bcVals[1]+(779.4228634059946*phiLx[0]-259.8076211353315*phiLy[0])*rdxCp2R3[0])*rdxCp2[1]+600.0*rdxCp2R4[0]*phiC[1]-692.8203230275509*rdxCp2R4[0]*bcVals[1])*omega-66150.0*phiC[1]*rdxCp2R4[1]-223020.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-600.0*rdxCp2R4[0]*phiC[1]))/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[2] = (((1153.545837840872*rdxCp2[0]*rdxCp2Sq[1]+2244.737846609264*rdxCp2Sq[0]*rdxCp2[1]+415.6921938165305*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-5455.960043841962*rho[0]*rdxCp2R3[1]-15328.64964698456*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-5715.767664977294*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-7482.459488697546*rdxCp2[0]*rdxCp2R3[1])-7482.459488697546*rdxCp2Sq[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(6547.152052610353*rdxCp2[0]*rdxCp2R3[1]+19381.64853669573*rdxCp2Sq[0]*rdxCp2Sq[1]+7205.331359486529*rdxCp2R3[0]*rdxCp2[1]+692.8203230275509*rdxCp2R4[0])*phiLx[3]+(32735.76026305177*rdxCp2R4[1]+113795.7380572752*rdxCp2[0]*rdxCp2R3[1]+64951.90528383289*rdxCp2Sq[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-18900.0*rdxCp2R4[1])-59400.0*rdxCp2[0]*rdxCp2R3[1]-28650.0*rdxCp2Sq[0]*rdxCp2Sq[1]-3000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(5670.0*rdxCp2[0]*rdxCp2R3[1]+16785.0*rdxCp2Sq[0]*rdxCp2Sq[1]+6240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiLx[2]+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[2]-16367.88013152588*phiLy[0]*rdxCp2R4[1]+((-6480.0*rdxCp2[0]*phiLy[1])-6300.0*rdxCp2[0]*phiLx[1]-3637.306695894642*rdxCp2[0]*bcVals[1]+((-51441.90898479563*phiLy[0])-5455.960043841962*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-6480.0*rdxCp2Sq[0]*phiLy[1])-8850.0*rdxCp2Sq[0]*phiLx[1]-20438.19952931275*rdxCp2Sq[0]*bcVals[1]+((-24811.62781842416*phiLy[0])-7664.324823492281*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-900.0*rdxCp2R3[0]*phiLy[1])-3000.0*rdxCp2R3[0]*phiLx[1]-7967.433714816835*rdxCp2R3[0]*bcVals[1]+((-2598.076211353316*phiLy[0])-2598.076211353316*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[2])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 
  phiC[3] = (((5670.0*rdxCp2R3[1]+15786.0*rdxCp2[0]*rdxCp2Sq[1]+4296.0*rdxCp2Sq[0]*rdxCp2[1]+240.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-4425.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-1800.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-18900.0*rdxCp2R4[1])-42120.0*rdxCp2[0]*rdxCp2R3[1]-11370.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-9450.0*rdxCp2[0]*rdxCp2R3[1])-25200.0*rdxCp2Sq[0]*rdxCp2Sq[1]-5000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-66150.0*rdxCp2R4[1])-223020.0*rdxCp2[0]*rdxCp2R3[1]-130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]-21240.0*rdxCp2R3[0]*rdxCp2[1]-600.0*rdxCp2R4[0])*phiC[3]+((-4950.0*rdxCp2[0]*rdxCp2R3[1])+10800.0*rdxCp2Sq[0]*rdxCp2Sq[1]+3000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-6235.382907247957*rdxCp2[0]*rdxCp2R3[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-8183.940065762942*rdxCp2[0]*rdxCp2R3[1])-21823.84017536785*rdxCp2Sq[0]*rdxCp2Sq[1]-4330.127018922193*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-16367.88013152588*phiLy[1]*rdxCp2R4[1]+((-36476.99000740053*rdxCp2[0]*phiLy[1])+9093.266739736604*rdxCp2[0]*phiLx[1]-10500.0*rdxCp2[0]*bcVals[1]+(7875.0*phiLx[0]-5400.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-9846.708841029065*rdxCp2Sq[0]*phiLy[1])-8850.0*rdxCp2Sq[0]*bcVals[1]-5400.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-519.6152422706631*rdxCp2R3[0]*phiLy[1])-866.0254037844386*rdxCp2R3[0]*phiLx[1]-2600.0*rdxCp2R3[0]*bcVals[1]+((-750.0*phiLy[0])-750.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0])*phiC[3])/(66150.0*rdxCp2R4[1]+223020.0*rdxCp2[0]*rdxCp2R3[1]+130335.0*rdxCp2Sq[0]*rdxCp2Sq[1]+21240.0*rdxCp2R3[0]*rdxCp2[1]+600.0*rdxCp2R4[0]); 

}

void MGpoissonDampedGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

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

  phiC[0] = ((216.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(51.96152422706631*rdxCp2Sq[1]+571.5767664977294*rdxCp2[0]*rdxCp2[1])*rho[2]+(571.5767664977294*rdxCp2[0]*rdxCp2[1]+51.96152422706631*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((300.0*rdxCp2[0]*rdxCp2Sq[1]+60.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(60.0*rdxCp2[0]*rdxCp2Sq[1]+300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(200.0*rdxCp2R3[1]+2220.0*rdxCp2[0]*rdxCp2Sq[1]+100.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(86.60254037844386*rdxCp2R3[1]+987.26896031426*rdxCp2[0]*rdxCp2Sq[1]+173.2050807568877*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(51.96152422706631*rdxCp2[0]*rdxCp2Sq[1]+259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(75.0*phiLy[0]-75.0*phiC[0])*rdxCp2R3[1]+(259.8076211353315*rdxCp2[0]*phiLy[1]+173.2050807568877*rdxCp2[0]*phiLx[1]+100.0*rdxCp2[0]*bcVals[1]+(855.0*phiLy[0]+150.0*phiLx[0]-1005.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2Sq[0]*phiLy[1]+987.26896031426*rdxCp2Sq[0]*phiLx[1]+2220.0*rdxCp2Sq[0]*bcVals[1]+(150.0*phiLy[0]+855.0*phiLx[0]-1005.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+86.60254037844386*rdxCp2R3[0]*phiLx[1]+200.0*rdxCp2R3[0]*bcVals[1]+(75.0*phiLx[0]-75.0*phiC[0])*rdxCp2R3[0])*omega+75.0*phiC[0]*rdxCp2R3[1]+1005.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+1005.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+75.0*phiC[0]*rdxCp2R3[0])/(75.0*rdxCp2R3[1]+1005.0*rdxCp2[0]*rdxCp2Sq[1]+1005.0*rdxCp2Sq[0]*rdxCp2[1]+75.0*rdxCp2R3[0]); 
  phiC[1] = (((155.8845726811989*rdxCp2Sq[1]+218.2384017536785*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0]*rho[0])*omega*volFac+((259.8076211353315*rdxCp2R3[1]+883.3459118601273*rdxCp2[0]*rdxCp2Sq[1]+103.9230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(259.8076211353315*rdxCp2Sq[0]*rdxCp2[1]-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(1991.858428704209*rdxCp2[0]*rdxCp2Sq[1]+86.60254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(750.0*rdxCp2[0]*rdxCp2Sq[1]+150.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(225.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(225.0*phiLy[1]-225.0*phiC[1])*rdxCp2R3[1]+(765.0*rdxCp2[0]*phiLy[1]-750.0*rdxCp2[0]*phiLx[1]-3015.0*rdxCp2[0]*phiC[1]+866.0254037844386*rdxCp2[0]*bcVals[1]+(649.5190528383289*phiLy[0]-649.5190528383289*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(90.0*rdxCp2Sq[0]*phiLy[1]-150.0*rdxCp2Sq[0]*phiLx[1]-3015.0*rdxCp2Sq[0]*phiC[1]+3031.088913245535*rdxCp2Sq[0]*bcVals[1]+(129.9038105676658*phiLy[0]-129.9038105676658*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-225.0*rdxCp2R3[0]*phiC[1]+259.8076211353315*rdxCp2R3[0]*bcVals[1])*omega+225.0*phiC[1]*rdxCp2R3[1]+3015.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+225.0*rdxCp2R3[0]*phiC[1])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[2] = (((218.2384017536785*rdxCp2[0]*rdxCp2[1]+155.8845726811989*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+129.9038105676658*rho[0]*rdxCp2Sq[1]+1428.941916244323*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((259.8076211353315*rdxCp2[0]*rdxCp2Sq[1]-259.8076211353315*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+883.3459118601273*rdxCp2Sq[0]*rdxCp2[1]+259.8076211353315*rdxCp2R3[0])*phiLx[3]+(259.8076211353315*rdxCp2R3[1]+3031.088913245535*rdxCp2[0]*rdxCp2Sq[1]+866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-150.0*rdxCp2[0]*rdxCp2Sq[1])-750.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(90.0*rdxCp2[0]*rdxCp2Sq[1]+765.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiLx[2]+((-225.0*rdxCp2R3[1])-3015.0*rdxCp2[0]*rdxCp2Sq[1]-3015.0*rdxCp2Sq[0]*rdxCp2[1]-225.0*rdxCp2R3[0])*phiC[2]+(225.0*rdxCp2[0]*phiLy[1]+150.0*rdxCp2[0]*phiLx[1]+86.60254037844386*rdxCp2[0]*bcVals[1]+(129.9038105676658*phiLx[0]-129.9038105676658*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-225.0*rdxCp2Sq[0]*phiLy[1])+750.0*rdxCp2Sq[0]*phiLx[1]+1991.858428704209*rdxCp2Sq[0]*bcVals[1]+(649.5190528383289*phiLx[0]-649.5190528383289*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0])*phiC[2])/(225.0*rdxCp2R3[1]+3015.0*rdxCp2[0]*rdxCp2Sq[1]+3015.0*rdxCp2Sq[0]*rdxCp2[1]+225.0*rdxCp2R3[0]); 
  phiC[3] = (((180.0*rdxCp2Sq[1]+1152.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*rho[3]+(363.7306695894642*rdxCp2[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0])*rho[2]+(259.8076211353315*rdxCp2Sq[1]+363.7306695894642*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-1500.0*rdxCp2[0]*rdxCp2Sq[1])-300.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-300.0*rdxCp2[0]*rdxCp2Sq[1])-1500.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-450.0*rdxCp2R3[1])-6030.0*rdxCp2[0]*rdxCp2Sq[1]-6030.0*rdxCp2Sq[0]*rdxCp2[1]-450.0*rdxCp2R3[0])*phiC[3]+(1300.0*rdxCp2[0]*rdxCp2Sq[1]+500.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(433.0127018922193*rdxCp2[0]*rdxCp2Sq[1]-433.0127018922193*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-259.8076211353315*rdxCp2[0]*rdxCp2Sq[1])-1299.038105676658*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-1299.038105676658*rdxCp2[0]*phiLy[1])-433.0127018922193*rdxCp2[0]*phiLx[1]+500.0*rdxCp2[0]*bcVals[1]+(375.0*phiLy[0]-375.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-259.8076211353315*rdxCp2Sq[0]*phiLy[1])+433.0127018922193*rdxCp2Sq[0]*phiLx[1]+1300.0*rdxCp2Sq[0]*bcVals[1]+(375.0*phiLx[0]-375.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0])*phiC[3])/(450.0*rdxCp2R3[1]+6030.0*rdxCp2[0]*rdxCp2Sq[1]+6030.0*rdxCp2Sq[0]*rdxCp2[1]+450.0*rdxCp2R3[0]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiLy[2]+(3.0*phiUy[0]+3.0*phiLy[0]-6.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-9.0*phiC[0]+12.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiLy[3]+(3.0*phiUy[1]+3.0*phiLy[1]-6.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-50.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]+8.660254037844386*phiC[0]-34.64101615137754*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiUx[2]-10.0*rdxCp2[1]*phiLy[2]+((-40.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[0]*phiUx[3]-10.0*rdxCp2[1]*phiLy[3]+((-40.0*rdxCp2[1])-50.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]+8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]-13.85640646055102*rdxCp2[0]*phiUx[1]-20.78460969082652*rdxCp2[0]*phiC[1]+(12.0*phiUx[0]-12.0*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.02777777777777778*(36.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiLy[3]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.08333333333333333*(12.0*rho[2]*volFac-13.85640646055102*rdxCp2[0]*phiUx[3]-20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+12.0*rdxCp2[0]*phiUx[2]-30.0*rdxCp2[1]*phiLy[2]+((-120.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+(25.98076211353316*phiUy[0]-25.98076211353316*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiUx[3]-30.0*rdxCp2[1]*phiLy[3]+((-120.0*rdxCp2[1])-60.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiLy[2]+(3.0*phiUy[0]+3.0*phiLy[0]-6.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]-3.464101615137754*rdxCp2[0]*phiC[1]+12.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiLy[3]+(3.0*phiUy[1]+3.0*phiLy[1]-6.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiLx[1]-50.0*rdxCp2[0]*phiC[1]+34.64101615137754*rdxCp2[0]*bcVals[1]+((-8.660254037844386*phiLx[0])-8.660254037844386*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac+3.464101615137754*rdxCp2[0]*phiLx[3]-3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]-10.0*rdxCp2[1]*phiLy[2]+3.0*rdxCp2[0]*phiLx[2]+((-40.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[1]*phiLy[3]-10.0*rdxCp2[0]*phiLx[3]+((-40.0*rdxCp2[1])-50.0*rdxCp2[0])*phiC[3]-8.660254037844386*rdxCp2[0]*phiLx[2]-8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]+13.85640646055102*rdxCp2[0]*phiLx[1]+20.78460969082652*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(12.0*phiLx[0]-12.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.02777777777777778*(36.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiLy[3]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.08333333333333333*(12.0*rho[2]*volFac+13.85640646055102*rdxCp2[0]*phiLx[3]+20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[1]*phiLy[2]+12.0*rdxCp2[0]*phiLx[2]+((-120.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+(25.98076211353316*phiUy[0]-25.98076211353316*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[1]*phiUy[3]-30.0*rdxCp2[1]*phiLy[3]-20.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-60.0*rdxCp2[0])*phiC[3]-17.32050807568877*rdxCp2[0]*phiLx[2]+17.32050807568877*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]-25.98076211353316*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiC[2]+12.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiLx[1]+(3.0*phiUx[0]+3.0*phiLx[0]-6.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-10.0*rdxCp2[0]*phiLx[1]-40.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiLx[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiUx[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-6.0*rdxCp2[0])*phiC[2]-34.64101615137754*rdxCp2[1]*bcVals[2]+(8.660254037844386*phiUy[0]+8.660254037844386*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[0]*phiUx[3]-10.0*rdxCp2[0]*phiLx[3]+((-50.0*rdxCp2[1])-40.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]-8.660254037844386*rdxCp2[0]*phiLx[2]+(8.660254037844386*phiUy[1]+8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-13.85640646055102*rdxCp2[1]*phiUy[2]-20.78460969082652*rdxCp2[1]*phiC[2]-8.0*rdxCp2[1]*bcVals[2]+(12.0*phiUy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.08333333333333333*(12.0*rho[1]*volFac-13.85640646055102*rdxCp2[1]*phiUy[3]-20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiUy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]-120.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.02777777777777778*(36.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiLx[3]-60.0*rdxCp2[1]*phiUy[2]+27.0*rdxCp2[0]*phiUx[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-30.0*rdxCp2[0]*phiUx[3]-30.0*rdxCp2[0]*phiLx[3]+((-60.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiLx[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac+12.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]-3.464101615137754*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiLx[1]+(3.0*phiUx[0]+3.0*phiLx[0]-6.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac+3.464101615137754*rdxCp2[1]*phiLy[3]-3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-10.0*rdxCp2[0]*phiLx[1]-40.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiLx[3]+34.64101615137754*rdxCp2[1]*bcVals[3]+3.0*rdxCp2[0]*phiUx[2]-10.0*rdxCp2[1]*phiLy[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-6.0*rdxCp2[0])*phiC[2]+((-8.660254037844386*phiLy[0])-8.660254037844386*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.25*(4.0*rho[3]*volFac-10.0*rdxCp2[0]*phiUx[3]-10.0*rdxCp2[1]*phiLy[3]-10.0*rdxCp2[0]*phiLx[3]+((-50.0*rdxCp2[1])-40.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]-8.660254037844386*rdxCp2[0]*phiLx[2]+((-8.660254037844386*phiLy[1])-8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac+8.0*rdxCp2[1]*bcVals[3]+13.85640646055102*rdxCp2[1]*phiLy[2]+20.78460969082652*rdxCp2[1]*phiC[2]+(12.0*phiLy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.08333333333333333*(12.0*rho[1]*volFac+13.85640646055102*rdxCp2[1]*phiLy[3]+20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiLy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]-120.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.02777777777777778*(36.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiLx[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+27.0*rdxCp2[0]*phiUx[2]-60.0*rdxCp2[1]*phiLy[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[0]*phiUx[3]-20.0*rdxCp2[1]*phiLy[3]-30.0*rdxCp2[0]*phiLx[3]+((-60.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]-25.98076211353316*rdxCp2[0]*phiLx[2]+(17.32050807568877*phiC[1]-17.32050807568877*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiC[2]+12.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-9.0*phiC[0]+12.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-50.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]+8.660254037844386*phiC[0]-34.64101615137754*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiUx[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]-34.64101615137754*rdxCp2[1]*bcVals[2]+(8.660254037844386*phiUy[0]+8.660254037844386*phiC[0])*rdxCp2[1]); 
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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-13.85640646055102*rdxCp2[1]*phiUy[2]-20.78460969082652*rdxCp2[1]*phiC[2]-8.0*rdxCp2[1]*bcVals[2]+(12.0*phiUy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-27.0*phiC[0]+36.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.08333333333333333*(12.0*rho[1]*volFac-13.85640646055102*rdxCp2[1]*phiUy[3]-20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiUy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]+25.98076211353316*phiC[0]-103.9230484541326*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.02777777777777778*(36.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiC[3]-60.0*rdxCp2[1]*phiUy[2]+27.0*rdxCp2[0]*phiUx[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-30.0*rdxCp2[0]*phiUx[3]+((-60.0*rdxCp2[1])-150.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiC[2]+36.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-27.0*phiC[0])*rdxCp2[1]-13.85640646055102*rdxCp2[0]*phiUx[1]-20.78460969082652*rdxCp2[0]*phiC[1]+(12.0*phiUx[0]-12.0*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.02777777777777778*(36.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiUy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.08333333333333333*(12.0*rho[2]*volFac-13.85640646055102*rdxCp2[0]*phiUx[3]-20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+12.0*rdxCp2[0]*phiUx[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]-103.9230484541326*rdxCp2[1]*bcVals[2]+(25.98076211353316*phiUy[0]+25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiUx[3]+((-150.0*rdxCp2[1])-60.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]+25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]-5.196152422706631*rdxCp2[1]*phiC[2]-2.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-3.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]-5.196152422706631*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-3.0*phiC[0]-2.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.05555555555555555*(18.0*rho[1]*volFac-20.78460969082652*rdxCp2[1]*phiUy[3]-31.17691453623978*rdxCp2[1]*phiC[3]+(18.0*phiUy[1]-18.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-90.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]-25.98076211353316*phiC[0]+34.64101615137754*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.05555555555555555*(18.0*rho[2]*volFac-20.78460969082652*rdxCp2[0]*phiUx[3]-31.17691453623978*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+18.0*rdxCp2[0]*phiUx[2]+((-90.0*rdxCp2[1])-18.0*rdxCp2[0])*phiC[2]+34.64101615137754*rdxCp2[1]*bcVals[2]+(25.98076211353316*phiUy[0]-25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.1666666666666667*(6.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[0]*phiUx[3]+((-30.0*rdxCp2[1])-30.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]-8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac+12.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]-3.464101615137754*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-9.0*phiC[0]+12.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac+3.464101615137754*rdxCp2[1]*phiLy[3]-3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-50.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]+8.660254037844386*phiC[0]-34.64101615137754*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiC[3]+34.64101615137754*rdxCp2[1]*bcVals[3]+3.0*rdxCp2[0]*phiUx[2]-10.0*rdxCp2[1]*phiLy[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+((-8.660254037844386*phiLy[0])-8.660254037844386*phiC[0])*rdxCp2[1]); 
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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac+8.0*rdxCp2[1]*bcVals[3]+13.85640646055102*rdxCp2[1]*phiLy[2]+20.78460969082652*rdxCp2[1]*phiC[2]+(12.0*phiLy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-27.0*phiC[0]+36.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.08333333333333333*(12.0*rho[1]*volFac+13.85640646055102*rdxCp2[1]*phiLy[3]+20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiLy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]+25.98076211353316*phiC[0]-103.9230484541326*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.02777777777777778*(36.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+27.0*rdxCp2[0]*phiUx[2]-60.0*rdxCp2[1]*phiLy[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[0]*phiUx[3]-20.0*rdxCp2[1]*phiLy[3]+((-60.0*rdxCp2[1])-150.0*rdxCp2[0])*phiC[3]+25.98076211353316*rdxCp2[0]*phiUx[2]+25.98076211353316*rdxCp2[0]*phiC[2]+(17.32050807568877*phiC[1]-17.32050807568877*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac+36.0*rdxCp2[1]*bcVals[3]+10.39230484541326*rdxCp2[1]*phiLy[2]-10.39230484541326*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-27.0*phiC[0])*rdxCp2[1]-13.85640646055102*rdxCp2[0]*phiUx[1]-20.78460969082652*rdxCp2[0]*phiC[1]+(12.0*phiUx[0]-12.0*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.02777777777777778*(36.0*rho[1]*volFac+31.17691453623978*rdxCp2[1]*phiLy[3]-31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiLy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.08333333333333333*(12.0*rho[2]*volFac-13.85640646055102*rdxCp2[0]*phiUx[3]-20.78460969082652*rdxCp2[0]*phiC[3]+103.9230484541326*rdxCp2[1]*bcVals[3]+12.0*rdxCp2[0]*phiUx[2]-30.0*rdxCp2[1]*phiLy[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+((-25.98076211353316*phiLy[0])-25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-20.0*rdxCp2[0]*phiUx[3]-30.0*rdxCp2[1]*phiLy[3]+((-150.0*rdxCp2[1])-60.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiC[2]+((-25.98076211353316*phiLy[1])-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac+2.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]+5.196152422706631*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-3.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]-5.196152422706631*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-3.0*phiC[0]-2.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.05555555555555555*(18.0*rho[1]*volFac+20.78460969082652*rdxCp2[1]*phiLy[3]+31.17691453623978*rdxCp2[1]*phiC[3]+(18.0*phiLy[1]-18.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-90.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]-25.98076211353316*phiC[0]+34.64101615137754*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.05555555555555555*(18.0*rho[2]*volFac-20.78460969082652*rdxCp2[0]*phiUx[3]-31.17691453623978*rdxCp2[0]*phiC[3]+34.64101615137754*rdxCp2[1]*bcVals[3]+18.0*rdxCp2[0]*phiUx[2]-30.0*rdxCp2[1]*phiLy[2]+((-90.0*rdxCp2[1])-18.0*rdxCp2[0])*phiC[2]+(25.98076211353316*phiC[0]-25.98076211353316*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.1666666666666667*(6.0*rho[3]*volFac-10.0*rdxCp2[0]*phiUx[3]-10.0*rdxCp2[1]*phiLy[3]+((-30.0*rdxCp2[1])-30.0*rdxCp2[0])*phiC[3]+8.660254037844386*rdxCp2[0]*phiUx[2]-8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiC[1]-8.660254037844386*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiC[2]+12.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]-3.464101615137754*rdxCp2[0]*phiC[1]+12.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiLx[1]-50.0*rdxCp2[0]*phiC[1]+34.64101615137754*rdxCp2[0]*bcVals[1]+((-8.660254037844386*phiLx[0])-8.660254037844386*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac+3.464101615137754*rdxCp2[0]*phiLx[3]-3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]-34.64101615137754*rdxCp2[1]*bcVals[2]+(8.660254037844386*phiUy[0]+8.660254037844386*phiC[0])*rdxCp2[1]); 
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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-13.85640646055102*rdxCp2[1]*phiUy[2]-20.78460969082652*rdxCp2[1]*phiC[2]-8.0*rdxCp2[1]*bcVals[2]+(12.0*phiUy[0]-12.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiLx[1]-10.39230484541326*rdxCp2[0]*phiC[1]+36.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-27.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.08333333333333333*(12.0*rho[1]*volFac-13.85640646055102*rdxCp2[1]*phiUy[3]-20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiUy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiLx[1]-150.0*rdxCp2[0]*phiC[1]+103.9230484541326*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiLx[0])-25.98076211353316*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.02777777777777778*(36.0*rho[2]*volFac+31.17691453623978*rdxCp2[0]*phiLx[3]-31.17691453623978*rdxCp2[0]*phiC[3]-60.0*rdxCp2[1]*phiUy[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-30.0*rdxCp2[0]*phiLx[3]+((-60.0*rdxCp2[1])-150.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiC[2]+36.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-27.0*phiC[0])*rdxCp2[1]+13.85640646055102*rdxCp2[0]*phiLx[1]+20.78460969082652*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(12.0*phiLx[0]-12.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.02777777777777778*(36.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiUy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.08333333333333333*(12.0*rho[2]*volFac+13.85640646055102*rdxCp2[0]*phiLx[3]+20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+12.0*rdxCp2[0]*phiLx[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]-103.9230484541326*rdxCp2[1]*bcVals[2]+(25.98076211353316*phiUy[0]+25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiLx[3]+((-150.0*rdxCp2[1])-60.0*rdxCp2[0])*phiC[3]-17.32050807568877*rdxCp2[0]*phiLx[2]+17.32050807568877*rdxCp2[0]*phiC[2]+(25.98076211353316*phiUy[1]+25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]-5.196152422706631*rdxCp2[1]*phiC[2]-2.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-3.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]+5.196152422706631*rdxCp2[0]*phiC[1]+2.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-3.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.05555555555555555*(18.0*rho[1]*volFac-20.78460969082652*rdxCp2[1]*phiUy[3]-31.17691453623978*rdxCp2[1]*phiC[3]+(18.0*phiUy[1]-18.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiLx[1]-90.0*rdxCp2[0]*phiC[1]+34.64101615137754*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiC[0]-25.98076211353316*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.05555555555555555*(18.0*rho[2]*volFac+20.78460969082652*rdxCp2[0]*phiLx[3]+31.17691453623978*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+18.0*rdxCp2[0]*phiLx[2]+((-90.0*rdxCp2[1])-18.0*rdxCp2[0])*phiC[2]+34.64101615137754*rdxCp2[1]*bcVals[2]+(25.98076211353316*phiUy[0]-25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.1666666666666667*(6.0*rho[3]*volFac-10.0*rdxCp2[1]*phiUy[3]-10.0*rdxCp2[0]*phiLx[3]+((-30.0*rdxCp2[1])-30.0*rdxCp2[0])*phiC[3]-8.660254037844386*rdxCp2[0]*phiLx[2]+8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiUy[1]-8.660254037844386*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.25*(4.0*rho[0]*volFac+12.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]-3.464101615137754*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]-3.464101615137754*rdxCp2[0]*phiC[1]+12.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.25*(4.0*rho[1]*volFac+3.464101615137754*rdxCp2[1]*phiLy[3]-3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiLx[1]-50.0*rdxCp2[0]*phiC[1]+34.64101615137754*rdxCp2[0]*bcVals[1]+((-8.660254037844386*phiLx[0])-8.660254037844386*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.25*(4.0*rho[2]*volFac+3.464101615137754*rdxCp2[0]*phiLx[3]-3.464101615137754*rdxCp2[0]*phiC[3]+34.64101615137754*rdxCp2[1]*bcVals[3]-10.0*rdxCp2[1]*phiLy[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+((-8.660254037844386*phiLy[0])-8.660254037844386*phiC[0])*rdxCp2[1]); 
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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac+8.0*rdxCp2[1]*bcVals[3]+13.85640646055102*rdxCp2[1]*phiLy[2]+20.78460969082652*rdxCp2[1]*phiC[2]+(12.0*phiLy[0]-12.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiLx[1]-10.39230484541326*rdxCp2[0]*phiC[1]+36.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-27.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.08333333333333333*(12.0*rho[1]*volFac+13.85640646055102*rdxCp2[1]*phiLy[3]+20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiLy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiLx[1]-150.0*rdxCp2[0]*phiC[1]+103.9230484541326*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiLx[0])-25.98076211353316*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.02777777777777778*(36.0*rho[2]*volFac+31.17691453623978*rdxCp2[0]*phiLx[3]-31.17691453623978*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]-60.0*rdxCp2[1]*phiLy[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-20.0*rdxCp2[1]*phiLy[3]-30.0*rdxCp2[0]*phiLx[3]+((-60.0*rdxCp2[1])-150.0*rdxCp2[0])*phiC[3]-25.98076211353316*rdxCp2[0]*phiLx[2]-25.98076211353316*rdxCp2[0]*phiC[2]+(17.32050807568877*phiC[1]-17.32050807568877*phiLy[1])*rdxCp2[1]); 

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

  resOut[0] = 0.08333333333333333*(12.0*rho[0]*volFac+36.0*rdxCp2[1]*bcVals[3]+10.39230484541326*rdxCp2[1]*phiLy[2]-10.39230484541326*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-27.0*phiC[0])*rdxCp2[1]+13.85640646055102*rdxCp2[0]*phiLx[1]+20.78460969082652*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(12.0*phiLx[0]-12.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.02777777777777778*(36.0*rho[1]*volFac+31.17691453623978*rdxCp2[1]*phiLy[3]-31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiLy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.08333333333333333*(12.0*rho[2]*volFac+13.85640646055102*rdxCp2[0]*phiLx[3]+20.78460969082652*rdxCp2[0]*phiC[3]+103.9230484541326*rdxCp2[1]*bcVals[3]-30.0*rdxCp2[1]*phiLy[2]+12.0*rdxCp2[0]*phiLx[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+((-25.98076211353316*phiLy[0])-25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.08333333333333333*(12.0*rho[3]*volFac-30.0*rdxCp2[1]*phiLy[3]-20.0*rdxCp2[0]*phiLx[3]+((-150.0*rdxCp2[1])-60.0*rdxCp2[0])*phiC[3]-17.32050807568877*rdxCp2[0]*phiLx[2]+17.32050807568877*rdxCp2[0]*phiC[2]+((-25.98076211353316*phiLy[1])-25.98076211353316*phiC[1])*rdxCp2[1]); 

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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac+2.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]+5.196152422706631*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-3.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]+5.196152422706631*rdxCp2[0]*phiC[1]+2.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-3.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.05555555555555555*(18.0*rho[1]*volFac+20.78460969082652*rdxCp2[1]*phiLy[3]+31.17691453623978*rdxCp2[1]*phiC[3]+(18.0*phiLy[1]-18.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiLx[1]-90.0*rdxCp2[0]*phiC[1]+34.64101615137754*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiC[0]-25.98076211353316*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.05555555555555555*(18.0*rho[2]*volFac+20.78460969082652*rdxCp2[0]*phiLx[3]+31.17691453623978*rdxCp2[0]*phiC[3]+34.64101615137754*rdxCp2[1]*bcVals[3]-30.0*rdxCp2[1]*phiLy[2]+18.0*rdxCp2[0]*phiLx[2]+((-90.0*rdxCp2[1])-18.0*rdxCp2[0])*phiC[2]+(25.98076211353316*phiC[0]-25.98076211353316*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.1666666666666667*(6.0*rho[3]*volFac-10.0*rdxCp2[1]*phiLy[3]-10.0*rdxCp2[0]*phiLx[3]+((-30.0*rdxCp2[1])-30.0*rdxCp2[0])*phiC[3]-8.660254037844386*rdxCp2[0]*phiLx[2]+8.660254037844386*rdxCp2[0]*phiC[2]+(8.660254037844386*phiC[1]-8.660254037844386*phiLy[1])*rdxCp2[1]); 

}

