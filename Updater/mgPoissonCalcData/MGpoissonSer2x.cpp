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
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = 1.0/dxLx[1]; 
  rdxUx[1]   = 1.0/dxUx[1]; 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = 1.0/dxLy[0]; 
  rdxUy[0]   = 1.0/dxUy[0]; 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = 1.0/dxLy[1]; 
  rdxUy[1]   = 1.0/dxUy[1]; 
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

  phiC[0] = ((((24600.0*rdxUx[0]-24600.0*rdxLx[0])*rdxUySq[1]+(24600.0*rdxUxSq[0]-24600.0*rdxLxSq[0])*rdxUy[1]+(24600.0*rdxLx[0]-24600.0*rdxUx[0])*rdxLySq[1]+(24600.0*rdxLxSq[0]-24600.0*rdxUxSq[0])*rdxLy[1])*rho[3]+((-1662.768775266122*rdxUyCu[1])+((-109742.739167564*rdxLy[1])-65332.95646149805*rdxUx[0]-65332.95646149805*rdxLx[0])*rdxUySq[1]+(109742.739167564*rdxLySq[1]-63670.18768623193*rdxUxSq[0]-19260.40498016591*rdxLx[0]*rdxUx[0]-63670.18768623193*rdxLxSq[0])*rdxUy[1]+1662.768775266122*rdxLyCu[1]+(65332.95646149805*rdxUx[0]+65332.95646149805*rdxLx[0])*rdxLySq[1]+(63670.18768623193*rdxUxSq[0]+19260.40498016591*rdxLx[0]*rdxUx[0]+63670.18768623193*rdxLxSq[0])*rdxLy[1])*rho[2]+((63670.18768623193*rdxLx[0]-63670.18768623193*rdxUx[0])*rdxUySq[1]+((19260.40498016591*rdxLx[0]-19260.40498016591*rdxUx[0])*rdxLy[1]-65332.95646149805*rdxUxSq[0]+65332.95646149805*rdxLxSq[0])*rdxUy[1]+(63670.18768623193*rdxLx[0]-63670.18768623193*rdxUx[0])*rdxLySq[1]+(65332.95646149805*rdxLxSq[0]-65332.95646149805*rdxUxSq[0])*rdxLy[1]-1662.768775266122*rdxUxCu[0]-109742.739167564*rdxLx[0]*rdxUxSq[0]+109742.739167564*rdxLxSq[0]*rdxUx[0]+1662.768775266122*rdxLxCu[0])*rho[1]+4416.0*rho[0]*rdxUyCu[1]+(300288.0*rho[0]*rdxLy[1]+(176968.0*rdxUx[0]+176968.0*rdxLx[0])*rho[0])*rdxUySq[1]+(300288.0*rho[0]*rdxLySq[1]+(578576.0*rdxUx[0]+578576.0*rdxLx[0])*rho[0]*rdxLy[1]+(176968.0*rdxUxSq[0]+578576.0*rdxLx[0]*rdxUx[0]+176968.0*rdxLxSq[0])*rho[0])*rdxUy[1]+4416.0*rho[0]*rdxLyCu[1]+(176968.0*rdxUx[0]+176968.0*rdxLx[0])*rho[0]*rdxLySq[1]+(176968.0*rdxUxSq[0]+578576.0*rdxLx[0]*rdxUx[0]+176968.0*rdxLxSq[0])*rho[0]*rdxLy[1]+(4416.0*rdxUxCu[0]+300288.0*rdxLx[0]*rdxUxSq[0]+300288.0*rdxLxSq[0]*rdxUx[0]+4416.0*rdxLxCu[0])*rho[0])*volFac+((47400.0*rdxUx[0]-47400.0*rdxLx[0])*rdxUyCu[1]+((20850.0*rdxUx[0]-20850.0*rdxLx[0])*rdxLy[1]+49200.0*rdxUxSq[0]-49200.0*rdxLxSq[0])*rdxUySq[1]+((90450.0*rdxUx[0]-90450.0*rdxLx[0])*rdxLySq[1]+(92250.0*rdxUxSq[0]-92250.0*rdxLxSq[0])*rdxLy[1]+1800.0*rdxUxCu[0]+118800.0*rdxLx[0]*rdxUxSq[0]-118800.0*rdxLxSq[0]*rdxUx[0]-1800.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(1800.0*rdxUx[0]*rdxUyCu[1]+(118800.0*rdxUx[0]*rdxLy[1]+49200.0*rdxUxSq[0]+92250.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-118800.0*rdxUx[0]*rdxLySq[1])+47400.0*rdxUxCu[0]+20850.0*rdxLx[0]*rdxUxSq[0]+90450.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-1800.0*rdxUx[0]*rdxLyCu[1]+((-49200.0*rdxUxSq[0])-92250.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-47400.0*rdxUxCu[0])-20850.0*rdxLx[0]*rdxUxSq[0]-90450.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((90450.0*rdxLx[0]-90450.0*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((20850.0*rdxLx[0]-20850.0*rdxUx[0])*rdxLySq[1]+(92250.0*rdxLxSq[0]-92250.0*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(47400.0*rdxLx[0]-47400.0*rdxUx[0])*rdxLyCu[1]+(49200.0*rdxLxSq[0]-49200.0*rdxUxSq[0])*rdxLySq[1]+((-1800.0*rdxUxCu[0])-118800.0*rdxLx[0]*rdxUxSq[0]+118800.0*rdxLxSq[0]*rdxUx[0]+1800.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-1800.0*rdxLx[0]*rdxUyCu[1])+((-118800.0*rdxLx[0]*rdxLy[1])-92250.0*rdxLx[0]*rdxUx[0]-49200.0*rdxLxSq[0])*rdxUySq[1]+(118800.0*rdxLx[0]*rdxLySq[1]-90450.0*rdxLx[0]*rdxUxSq[0]-20850.0*rdxLxSq[0]*rdxUx[0]-47400.0*rdxLxCu[0])*rdxUy[1]+1800.0*rdxLx[0]*rdxLyCu[1]+(92250.0*rdxLx[0]*rdxUx[0]+49200.0*rdxLxSq[0])*rdxLySq[1]+(90450.0*rdxLx[0]*rdxUxSq[0]+20850.0*rdxLxSq[0]*rdxUx[0]+47400.0*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((-1662.768775266122*rdxUyR4[1])+((-114523.1993964541*rdxLy[1])-67203.57133367243*rdxUx[0]-67203.57133367243*rdxLx[0])*rdxUyCu[1]+((-210548.0961680727*rdxLySq[1])+((-313163.4462624908*rdxUx[0])-313163.4462624908*rdxLx[0])*rdxLy[1]-67931.03267285136*rdxUxSq[0]-304737.0190836682*rdxLx[0]*rdxUx[0]-67931.03267285136*rdxLxSq[0])*rdxUySq[1]+((-3117.691453623978*rdxLyCu[1])+((-124369.9082374832*rdxUx[0])-124369.9082374832*rdxLx[0])*rdxLySq[1]+((-123642.4468983043*rdxUxSq[0])-321589.8734413133*rdxLx[0]*rdxUx[0]-123642.4468983043*rdxLxSq[0])*rdxLy[1]-2390.23011444505*rdxUxCu[0]-162535.6477822634*rdxLx[0]*rdxUxSq[0]-162535.6477822634*rdxLxSq[0]*rdxUx[0]-2390.23011444505*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-1870.614872174387*rdxUx[0]*rdxUyCu[1])+((-123460.5815635095*rdxUx[0]*rdxLy[1])-46869.29485281381*rdxUxSq[0]-100129.8571855568*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(123460.5815635095*rdxUx[0]*rdxLySq[1]-44998.67998063943*rdxUxCu[0]-21667.95560268665*rdxLx[0]*rdxUxSq[0]-98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+1870.614872174387*rdxUx[0]*rdxLyCu[1]+(46869.29485281381*rdxUxSq[0]+100129.8571855568*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(44998.67998063943*rdxUxCu[0]+21667.95560268665*rdxLx[0]*rdxUxSq[0]+98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+(3117.691453623978*rdxLy[1]*rdxUyCu[1]+(210548.0961680727*rdxLySq[1]+(124369.9082374832*rdxUx[0]+124369.9082374832*rdxLx[0])*rdxLy[1])*rdxUySq[1]+(114523.1993964541*rdxLyCu[1]+(313163.4462624908*rdxUx[0]+313163.4462624908*rdxLx[0])*rdxLySq[1]+(123642.4468983043*rdxUxSq[0]+321589.8734413133*rdxLx[0]*rdxUx[0]+123642.4468983043*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+1662.768775266122*rdxLyR4[1]+(67203.57133367243*rdxUx[0]+67203.57133367243*rdxLx[0])*rdxLyCu[1]+(67931.03267285136*rdxUxSq[0]+304737.0190836682*rdxLx[0]*rdxUx[0]+67931.03267285136*rdxLxSq[0])*rdxLySq[1]+(2390.23011444505*rdxUxCu[0]+162535.6477822634*rdxLx[0]*rdxUxSq[0]+162535.6477822634*rdxLxSq[0]*rdxUx[0]+2390.23011444505*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-1870.614872174387*rdxLx[0]*rdxUyCu[1])+((-123460.5815635095*rdxLx[0]*rdxLy[1])-100129.8571855568*rdxLx[0]*rdxUx[0]-46869.29485281381*rdxLxSq[0])*rdxUySq[1]+(123460.5815635095*rdxLx[0]*rdxLySq[1]-98259.24231338239*rdxLx[0]*rdxUxSq[0]-21667.95560268665*rdxLxSq[0]*rdxUx[0]-44998.67998063943*rdxLxCu[0])*rdxUy[1]+1870.614872174387*rdxLx[0]*rdxLyCu[1]+(100129.8571855568*rdxLx[0]*rdxUx[0]+46869.29485281381*rdxLxSq[0])*rdxLySq[1]+(98259.24231338239*rdxLx[0]*rdxUxSq[0]+21667.95560268665*rdxLxSq[0]*rdxUx[0]+44998.67998063943*rdxLxCu[0])*rdxLy[1])*phiLx[2]+1584.0*phiUy[0]*rdxUyR4[1]+((109512.0*phiUy[0]+3384.0*phiLy[0])*rdxLy[1]+(44998.67998063943*rdxLx[0]-44998.67998063943*rdxUx[0])*phiUy[1]-2390.23011444505*rdxUx[0]*phiUx[1]+2390.23011444505*rdxLx[0]*phiLx[1]+(64182.0*phiUy[0]+2484.0*phiUx[0])*rdxUx[0]+(64182.0*phiUy[0]+2484.0*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((228312.0*phiUy[0]+228312.0*phiLy[0])*rdxLySq[1]+((21667.95560268665*rdxLx[0]-21667.95560268665*rdxUx[0])*phiUy[1]-162535.6477822634*rdxUx[0]*phiUx[1]+(98259.24231338239*rdxLx[0]-98259.24231338239*rdxUx[0])*phiLy[1]+162535.6477822634*rdxLx[0]*phiLx[1]+(325449.0*phiUy[0]+168912.0*phiUx[0]+134907.0*phiLy[0])*rdxUx[0]+(325449.0*phiUy[0]+134907.0*phiLy[0]+168912.0*phiLx[0])*rdxLx[0])*rdxLy[1]+(46869.29485281381*rdxLxSq[0]-46869.29485281381*rdxUxSq[0])*phiUy[1]+((-67931.03267285136*rdxUxSq[0])-123642.4468983043*rdxLx[0]*rdxUx[0])*phiUx[1]+(123642.4468983043*rdxLx[0]*rdxUx[0]+67931.03267285136*rdxLxSq[0])*phiLx[1]+(65082.0*phiUy[0]+65082.0*phiUx[0])*rdxUxSq[0]+(315024.0*phiUy[0]+134007.0*phiUx[0]+134007.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(65082.0*phiUy[0]+65082.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((3384.0*phiUy[0]+109512.0*phiLy[0])*rdxLyCu[1]+((98259.24231338239*rdxLx[0]-98259.24231338239*rdxUx[0])*phiUy[1]-162535.6477822634*rdxUx[0]*phiUx[1]+(21667.95560268665*rdxLx[0]-21667.95560268665*rdxUx[0])*phiLy[1]+162535.6477822634*rdxLx[0]*phiLx[1]+(134907.0*phiUy[0]+168912.0*phiUx[0]+325449.0*phiLy[0])*rdxUx[0]+(134907.0*phiUy[0]+325449.0*phiLy[0]+168912.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((100129.8571855568*rdxLxSq[0]-100129.8571855568*rdxUxSq[0])*phiUy[1]+((-304737.0190836682*rdxUxSq[0])-321589.8734413133*rdxLx[0]*rdxUx[0])*phiUx[1]+(100129.8571855568*rdxLxSq[0]-100129.8571855568*rdxUxSq[0])*phiLy[1]+(321589.8734413133*rdxLx[0]*rdxUx[0]+304737.0190836682*rdxLxSq[0])*phiLx[1]+(134007.0*phiUy[0]+315024.0*phiUx[0]+134007.0*phiLy[0])*rdxUxSq[0]+(335874.0*phiUy[0]+335874.0*phiUx[0]+335874.0*phiLy[0]+335874.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(134007.0*phiUy[0]+134007.0*phiLy[0]+315024.0*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-1870.614872174387*rdxUxCu[0])-123460.5815635095*rdxLx[0]*rdxUxSq[0]+123460.5815635095*rdxLxSq[0]*rdxUx[0]+1870.614872174387*rdxLxCu[0])*phiUy[1]+((-67203.57133367243*rdxUxCu[0])-313163.4462624908*rdxLx[0]*rdxUxSq[0]-124369.9082374832*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(124369.9082374832*rdxLx[0]*rdxUxSq[0]+313163.4462624908*rdxLxSq[0]*rdxUx[0]+67203.57133367243*rdxLxCu[0])*phiLx[1]+(2484.0*phiUy[0]+64182.0*phiUx[0])*rdxUxCu[0]+(168912.0*phiUy[0]+325449.0*phiUx[0]+134907.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(168912.0*phiUy[0]+134907.0*phiUx[0]+325449.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(2484.0*phiUy[0]+64182.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+1584.0*phiLy[0]*rdxLyR4[1]+((-2390.23011444505*rdxUx[0]*phiUx[1])+(44998.67998063943*rdxLx[0]-44998.67998063943*rdxUx[0])*phiLy[1]+2390.23011444505*rdxLx[0]*phiLx[1]+(2484.0*phiUx[0]+64182.0*phiLy[0])*rdxUx[0]+(64182.0*phiLy[0]+2484.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+(((-67931.03267285136*rdxUxSq[0])-123642.4468983043*rdxLx[0]*rdxUx[0])*phiUx[1]+(46869.29485281381*rdxLxSq[0]-46869.29485281381*rdxUxSq[0])*phiLy[1]+(123642.4468983043*rdxLx[0]*rdxUx[0]+67931.03267285136*rdxLxSq[0])*phiLx[1]+(65082.0*phiUx[0]+65082.0*phiLy[0])*rdxUxSq[0]+(134007.0*phiUx[0]+315024.0*phiLy[0]+134007.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(65082.0*phiLy[0]+65082.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-67203.57133367243*rdxUxCu[0])-313163.4462624908*rdxLx[0]*rdxUxSq[0]-124369.9082374832*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-1870.614872174387*rdxUxCu[0])-123460.5815635095*rdxLx[0]*rdxUxSq[0]+123460.5815635095*rdxLxSq[0]*rdxUx[0]+1870.614872174387*rdxLxCu[0])*phiLy[1]+(124369.9082374832*rdxLx[0]*rdxUxSq[0]+313163.4462624908*rdxLxSq[0]*rdxUx[0]+67203.57133367243*rdxLxCu[0])*phiLx[1]+(64182.0*phiUx[0]+2484.0*phiLy[0])*rdxUxCu[0]+(325449.0*phiUx[0]+168912.0*phiLy[0]+134907.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(134907.0*phiUx[0]+168912.0*phiLy[0]+325449.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(2484.0*phiLy[0]+64182.0*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-1662.768775266122*rdxUxR4[0])-114523.1993964541*rdxLx[0]*rdxUxCu[0]-210548.0961680727*rdxLxSq[0]*rdxUxSq[0]-3117.691453623978*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(3117.691453623978*rdxLx[0]*rdxUxCu[0]+210548.0961680727*rdxLxSq[0]*rdxUxSq[0]+114523.1993964541*rdxLxCu[0]*rdxUx[0]+1662.768775266122*rdxLxR4[0])*phiLx[1]+1584.0*phiUx[0]*rdxUxR4[0]+(109512.0*phiUx[0]+3384.0*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(228312.0*phiUx[0]+228312.0*phiLx[0])*rdxLxSq[0]*rdxUxSq[0]+(3384.0*phiUx[0]+109512.0*phiLx[0])*rdxLxCu[0]*rdxUx[0]+1584.0*phiLx[0]*rdxLxR4[0])/(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0]); 
  phiC[1] = -(1.0*(((831.384387633061*rdxUyCu[1]+(54871.36958378201*rdxLy[1]+25565.06991971662*rdxUx[0]+25565.06991971662*rdxLx[0])*rdxUySq[1]+((-54871.36958378201*rdxLySq[1])+24733.68553208356*rdxUxSq[0]-4572.614131981835*rdxLx[0]*rdxUx[0]+24733.68553208356*rdxLxSq[0])*rdxUy[1]-831.384387633061*rdxLyCu[1]+((-25565.06991971662*rdxUx[0])-25565.06991971662*rdxLx[0])*rdxLySq[1]+((-24733.68553208356*rdxUxSq[0])+4572.614131981835*rdxLx[0]*rdxUx[0]-24733.68553208356*rdxLxSq[0])*rdxLy[1])*rho[3]+((63960.0*rdxLx[0]-63960.0*rdxUx[0])*rdxUySq[1]+(63960.0*rdxLxSq[0]-63960.0*rdxUxSq[0])*rdxUy[1]+(63960.0*rdxUx[0]-63960.0*rdxLx[0])*rdxLySq[1]+(63960.0*rdxUxSq[0]-63960.0*rdxLxSq[0])*rdxLy[1])*rho[2]+((-2208.0*rdxUyCu[1])+((-150144.0*rdxLy[1])-70104.0*rdxUx[0]-70104.0*rdxLx[0])*rdxUySq[1]+((-150144.0*rdxLySq[1])+((-283728.0*rdxUx[0])-283728.0*rdxLx[0])*rdxLy[1]-69624.0*rdxUxSq[0]-251568.0*rdxLx[0]*rdxUx[0]-69624.0*rdxLxSq[0])*rdxUy[1]-2208.0*rdxLyCu[1]+((-70104.0*rdxUx[0])-70104.0*rdxLx[0])*rdxLySq[1]+((-69624.0*rdxUxSq[0])-251568.0*rdxLx[0]*rdxUx[0]-69624.0*rdxLxSq[0])*rdxLy[1]-1728.0*rdxUxCu[0]-117504.0*rdxLx[0]*rdxUxSq[0]-117504.0*rdxLxSq[0]*rdxUx[0]-1728.0*rdxLxCu[0])*rho[1]+(165542.487984203*rdxUx[0]-165542.487984203*rdxLx[0])*rho[0]*rdxUySq[1]+((50077.05294843137*rdxUx[0]-50077.05294843137*rdxLx[0])*rho[0]*rdxLy[1]+(169865.6867998949*rdxUxSq[0]-169865.6867998949*rdxLxSq[0])*rho[0])*rdxUy[1]+(165542.487984203*rdxUx[0]-165542.487984203*rdxLx[0])*rho[0]*rdxLySq[1]+(169865.6867998949*rdxUxSq[0]-169865.6867998949*rdxLxSq[0])*rho[0]*rdxLy[1]+(4323.198815691917*rdxUxCu[0]+285331.1218356665*rdxLx[0]*rdxUxSq[0]-285331.1218356665*rdxLxSq[0]*rdxUx[0]-4323.198815691917*rdxLxCu[0])*rho[0])*volFac+(1662.768775266122*rdxUyR4[1]+(114523.1993964541*rdxLy[1]+53520.3699538783*rdxUx[0]+53520.3699538783*rdxLx[0])*rdxUyCu[1]+(210548.0961680727*rdxLySq[1]+(307144.569706189*rdxUx[0]+307144.569706189*rdxLx[0])*rdxLy[1]+53728.21605078656*rdxUxSq[0]+276331.3858395386*rdxLx[0]*rdxUx[0]+53728.21605078656*rdxLxSq[0])*rdxUySq[1]+(3117.691453623978*rdxLyCu[1]+(98259.24231338239*rdxUx[0]+98259.24231338239*rdxLx[0])*rdxLySq[1]+(97012.16573193281*rdxUxSq[0]+268329.3111085704*rdxLx[0]*rdxUx[0]+97012.16573193281*rdxLxSq[0])*rdxLy[1]+1870.614872174387*rdxUxCu[0]+127201.8113078583*rdxLx[0]*rdxUxSq[0]+127201.8113078583*rdxLxSq[0]*rdxUx[0]+1870.614872174387*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-727.4613391789284*rdxUx[0]*rdxUyCu[1])+((-48012.44838580926*rdxUx[0]*rdxLy[1])+46869.29485281381*rdxUxSq[0]-91608.1672123179*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(48012.44838580926*rdxUx[0]*rdxLySq[1]+47596.75619199274*rdxUxCu[0]+4001.037365484106*rdxLx[0]*rdxUxSq[0]-90880.70587313897*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+727.4613391789284*rdxUx[0]*rdxLyCu[1]+(91608.1672123179*rdxLx[0]*rdxUx[0]-46869.29485281381*rdxUxSq[0])*rdxLySq[1]+((-47596.75619199274*rdxUxCu[0])-4001.037365484106*rdxLx[0]*rdxUxSq[0]+90880.70587313897*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((-3117.691453623978*rdxLy[1]*rdxUyCu[1])+(((-98259.24231338239*rdxUx[0])-98259.24231338239*rdxLx[0])*rdxLy[1]-210548.0961680727*rdxLySq[1])*rdxUySq[1]+((-114523.1993964541*rdxLyCu[1])+((-307144.569706189*rdxUx[0])-307144.569706189*rdxLx[0])*rdxLySq[1]+((-97012.16573193281*rdxUxSq[0])-268329.3111085704*rdxLx[0]*rdxUx[0]-97012.16573193281*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-1662.768775266122*rdxLyR4[1]+((-53520.3699538783*rdxUx[0])-53520.3699538783*rdxLx[0])*rdxLyCu[1]+((-53728.21605078656*rdxUxSq[0])-276331.3858395386*rdxLx[0]*rdxUx[0]-53728.21605078656*rdxLxSq[0])*rdxLySq[1]+((-1870.614872174387*rdxUxCu[0])-127201.8113078583*rdxLx[0]*rdxUxSq[0]-127201.8113078583*rdxLxSq[0]*rdxUx[0]-1870.614872174387*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-727.4613391789284*rdxLx[0]*rdxUyCu[1])+((-48012.44838580926*rdxLx[0]*rdxLy[1])-91608.1672123179*rdxLx[0]*rdxUx[0]+46869.29485281381*rdxLxSq[0])*rdxUySq[1]+(48012.44838580926*rdxLx[0]*rdxLySq[1]-90880.70587313897*rdxLx[0]*rdxUxSq[0]+4001.037365484106*rdxLxSq[0]*rdxUx[0]+47596.75619199274*rdxLxCu[0])*rdxUy[1]+727.4613391789284*rdxLx[0]*rdxLyCu[1]+(91608.1672123179*rdxLx[0]*rdxUx[0]-46869.29485281381*rdxLxSq[0])*rdxLySq[1]+(90880.70587313897*rdxLx[0]*rdxUxSq[0]-4001.037365484106*rdxLxSq[0]*rdxUx[0]-47596.75619199274*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((61620.0*rdxLx[0]-61620.0*rdxUx[0])*rdxUyCu[1]+((27105.0*rdxLx[0]-27105.0*rdxUx[0])*rdxLy[1]-63960.0*rdxUxSq[0]+63960.0*rdxLxSq[0])*rdxUySq[1]+((117585.0*rdxLx[0]-117585.0*rdxUx[0])*rdxLySq[1]+(119925.0*rdxLxSq[0]-119925.0*rdxUxSq[0])*rdxLy[1]-2340.0*rdxUxCu[0]-154440.0*rdxLx[0]*rdxUxSq[0]+154440.0*rdxLxSq[0]*rdxUx[0]+2340.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(900.0*rdxUx[0]*rdxUyCu[1]+(59400.0*rdxUx[0]*rdxLy[1]-44280.0*rdxUxSq[0]+99630.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-59400.0*rdxUx[0]*rdxLySq[1])-45180.0*rdxUxCu[0]-4950.0*rdxLx[0]*rdxUxSq[0]+98730.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-900.0*rdxUx[0]*rdxLyCu[1]+(44280.0*rdxUxSq[0]-99630.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(45180.0*rdxUxCu[0]+4950.0*rdxLx[0]*rdxUxSq[0]-98730.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+((117585.0*rdxUx[0]-117585.0*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((27105.0*rdxUx[0]-27105.0*rdxLx[0])*rdxLySq[1]+(119925.0*rdxUxSq[0]-119925.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(61620.0*rdxUx[0]-61620.0*rdxLx[0])*rdxLyCu[1]+(63960.0*rdxUxSq[0]-63960.0*rdxLxSq[0])*rdxLySq[1]+(2340.0*rdxUxCu[0]+154440.0*rdxLx[0]*rdxUxSq[0]-154440.0*rdxLxSq[0]*rdxUx[0]-2340.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-900.0*rdxLx[0]*rdxUyCu[1])+((-59400.0*rdxLx[0]*rdxLy[1])-99630.0*rdxLx[0]*rdxUx[0]+44280.0*rdxLxSq[0])*rdxUySq[1]+(59400.0*rdxLx[0]*rdxLySq[1]-98730.0*rdxLx[0]*rdxUxSq[0]+4950.0*rdxLxSq[0]*rdxUx[0]+45180.0*rdxLxCu[0])*rdxUy[1]+900.0*rdxLx[0]*rdxLyCu[1]+(99630.0*rdxLx[0]*rdxUx[0]-44280.0*rdxLxSq[0])*rdxLySq[1]+(98730.0*rdxLx[0]*rdxUxSq[0]-4950.0*rdxLxSq[0]*rdxUx[0]-45180.0*rdxLxCu[0])*rdxLy[1])*phiLx[2]-1584.0*phiUy[1]*rdxUyR4[1]+(((-109512.0*phiUy[1])-3384.0*phiLy[1])*rdxLy[1]+((-51192.0*rdxUx[0])-51192.0*rdxLx[0])*phiUy[1]+966.0*rdxUx[0]*phiUx[1]+966.0*rdxLx[0]*phiLx[1]+(58498.28397483125*phiUy[0]-1195.115057222525*phiUx[0])*rdxUx[0]+(1195.115057222525*phiLx[0]-58498.28397483125*phiUy[0])*rdxLx[0])*rdxUyCu[1]+(((-228312.0*phiUy[1])-228312.0*phiLy[1])*rdxLySq[1]+(((-319194.0*rdxUx[0])-319194.0*rdxLx[0])*phiUy[1]+65688.0*rdxUx[0]*phiUx[1]+((-106542.0*rdxUx[0])-106542.0*rdxLx[0])*phiLy[1]+65688.0*rdxLx[0]*phiLx[1]+(28168.34228349264*phiUy[0]-81267.82389113172*phiUx[0]+127737.0150073971*phiLy[0])*rdxUx[0]+((-28168.34228349264*phiUy[0])-127737.0150073971*phiLy[0]+81267.82389113172*phiLx[0])*rdxLx[0])*rdxLy[1]+((-51552.0*rdxUxSq[0])-287964.0*rdxLx[0]*rdxUx[0]-51552.0*rdxLxSq[0])*phiUy[1]+(120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxUxSq[0])*phiUx[1]+(120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*phiLx[1]+(60930.08330865796*phiUy[0]+55172.74642429901*phiUx[0])*rdxUxSq[0]+(131062.5525579294*phiLx[0]-131062.5525579294*phiUx[0])*rdxLx[0]*rdxUx[0]+((-60930.08330865796*phiUy[0])-55172.74642429901*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-3384.0*phiUy[1])-109512.0*phiLy[1])*rdxLyCu[1]+(((-106542.0*rdxUx[0])-106542.0*rdxLx[0])*phiUy[1]+65688.0*rdxUx[0]*phiUx[1]+((-319194.0*rdxUx[0])-319194.0*rdxLx[0])*phiLy[1]+65688.0*rdxLx[0]*phiLx[1]+(127737.0150073971*phiUy[0]-81267.82389113172*phiUx[0]+28168.34228349264*phiLy[0])*rdxUx[0]+((-127737.0150073971*phiUy[0])-28168.34228349264*phiLy[0]+81267.82389113172*phiLx[0])*rdxLx[0])*rdxLySq[1]+(((-105102.0*rdxUxSq[0])-278064.0*rdxLx[0]*rdxUx[0]-105102.0*rdxLxSq[0])*phiUy[1]+(97026.0*rdxUxSq[0]+151236.0*rdxLx[0]*rdxUx[0])*phiUx[1]+((-105102.0*rdxUxSq[0])-278064.0*rdxLx[0]*rdxUx[0]-105102.0*rdxLxSq[0])*phiLy[1]+(151236.0*rdxLx[0]*rdxUx[0]+97026.0*rdxLxSq[0])*phiLx[1]+(130168.8143412238*phiUy[0]-125403.9425696018*phiUx[0]+130168.8143412238*phiLy[0])*rdxUxSq[0]+(181740.6271365871*phiLx[0]-181740.6271365871*phiUx[0])*rdxLx[0]*rdxUx[0]+((-130168.8143412238*phiUy[0])-130168.8143412238*phiLy[0]+125403.9425696018*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-1944.0*rdxUxCu[0])-132192.0*rdxLx[0]*rdxUxSq[0]-132192.0*rdxLxSq[0]*rdxUx[0]-1944.0*rdxLxCu[0])*phiUy[1]+((-61482.0*rdxUxCu[0])+110061.0*rdxLx[0]*rdxUxSq[0]+122403.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(122403.0*rdxLx[0]*rdxUxSq[0]+110061.0*rdxLxSq[0]*rdxUx[0]-61482.0*rdxLxCu[0])*phiLx[1]+(2431.799333826703*phiUy[0]+57864.35337926103*phiUx[0])*rdxUxCu[0]+(160498.7560325623*phiUy[0]-136165.1742370272*phiUx[0]+133234.5442706207*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-160498.7560325623*phiUy[0])-133234.5442706207*phiUx[0]+136165.1742370272*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-2431.799333826703*phiUy[0])-57864.35337926103*phiLx[0])*rdxLxCu[0])*rdxUy[1]-1584.0*phiLy[1]*rdxLyR4[1]+(966.0*rdxUx[0]*phiUx[1]+((-51192.0*rdxUx[0])-51192.0*rdxLx[0])*phiLy[1]+966.0*rdxLx[0]*phiLx[1]+(58498.28397483125*phiLy[0]-1195.115057222525*phiUx[0])*rdxUx[0]+(1195.115057222525*phiLx[0]-58498.28397483125*phiLy[0])*rdxLx[0])*rdxLyCu[1]+((120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxUxSq[0])*phiUx[1]+((-51552.0*rdxUxSq[0])-287964.0*rdxLx[0]*rdxUx[0]-51552.0*rdxLxSq[0])*phiLy[1]+(120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*phiLx[1]+(55172.74642429901*phiUx[0]+60930.08330865796*phiLy[0])*rdxUxSq[0]+(131062.5525579294*phiLx[0]-131062.5525579294*phiUx[0])*rdxLx[0]*rdxUx[0]+((-60930.08330865796*phiLy[0])-55172.74642429901*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-61482.0*rdxUxCu[0])+110061.0*rdxLx[0]*rdxUxSq[0]+122403.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-1944.0*rdxUxCu[0])-132192.0*rdxLx[0]*rdxUxSq[0]-132192.0*rdxLxSq[0]*rdxUx[0]-1944.0*rdxLxCu[0])*phiLy[1]+(122403.0*rdxLx[0]*rdxUxSq[0]+110061.0*rdxLxSq[0]*rdxUx[0]-61482.0*rdxLxCu[0])*phiLx[1]+(57864.35337926103*phiUx[0]+2431.799333826703*phiLy[0])*rdxUxCu[0]+((-136165.1742370272*phiUx[0])+160498.7560325623*phiLy[0]+133234.5442706207*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-133234.5442706207*phiUx[0])-160498.7560325623*phiLy[0]+136165.1742370272*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-2431.799333826703*phiLy[0])-57864.35337926103*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-1584.0*rdxUxR4[0])-103032.0*rdxLx[0]*rdxUxCu[0]+205848.0*rdxLxSq[0]*rdxUxSq[0]+3096.0*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(3096.0*rdxLx[0]*rdxUxCu[0]+205848.0*rdxLxSq[0]*rdxUxSq[0]-103032.0*rdxLxCu[0]*rdxUx[0]-1584.0*rdxLxR4[0])*phiLx[1]+1496.491897739509*phiUx[0]*rdxUxR4[0]+(96897.85037863323*phiUx[0]+3367.106769913895*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(224099.6616864915*phiLx[0]-224099.6616864915*phiUx[0])*rdxLxSq[0]*rdxUxSq[0]+((-3367.106769913895*phiUx[0])-96897.85037863323*phiLx[0])*rdxLxCu[0]*rdxUx[0]-1496.491897739509*phiLx[0]*rdxLxR4[0]))/(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0]); 
  phiC[2] = -(1.0*((((24733.68553208356*rdxUx[0]-24733.68553208356*rdxLx[0])*rdxUySq[1]+((4572.614131981835*rdxLx[0]-4572.614131981835*rdxUx[0])*rdxLy[1]+25565.06991971662*rdxUxSq[0]-25565.06991971662*rdxLxSq[0])*rdxUy[1]+(24733.68553208356*rdxUx[0]-24733.68553208356*rdxLx[0])*rdxLySq[1]+(25565.06991971662*rdxUxSq[0]-25565.06991971662*rdxLxSq[0])*rdxLy[1]+831.384387633061*rdxUxCu[0]+54871.36958378201*rdxLx[0]*rdxUxSq[0]-54871.36958378201*rdxLxSq[0]*rdxUx[0]-831.384387633061*rdxLxCu[0])*rho[3]+((-1728.0*rdxUyCu[1])+((-117504.0*rdxLy[1])-69624.0*rdxUx[0]-69624.0*rdxLx[0])*rdxUySq[1]+((-117504.0*rdxLySq[1])+((-251568.0*rdxUx[0])-251568.0*rdxLx[0])*rdxLy[1]-70104.0*rdxUxSq[0]-283728.0*rdxLx[0]*rdxUx[0]-70104.0*rdxLxSq[0])*rdxUy[1]-1728.0*rdxLyCu[1]+((-69624.0*rdxUx[0])-69624.0*rdxLx[0])*rdxLySq[1]+((-70104.0*rdxUxSq[0])-283728.0*rdxLx[0]*rdxUx[0]-70104.0*rdxLxSq[0])*rdxLy[1]-2208.0*rdxUxCu[0]-150144.0*rdxLx[0]*rdxUxSq[0]-150144.0*rdxLxSq[0]*rdxUx[0]-2208.0*rdxLxCu[0])*rho[2]+((63960.0*rdxLx[0]-63960.0*rdxUx[0])*rdxUySq[1]+(63960.0*rdxLxSq[0]-63960.0*rdxUxSq[0])*rdxUy[1]+(63960.0*rdxUx[0]-63960.0*rdxLx[0])*rdxLySq[1]+(63960.0*rdxUxSq[0]-63960.0*rdxLxSq[0])*rdxLy[1])*rho[1]+4323.198815691917*rho[0]*rdxUyCu[1]+(285331.1218356665*rho[0]*rdxLy[1]+(169865.6867998949*rdxUx[0]+169865.6867998949*rdxLx[0])*rho[0])*rdxUySq[1]+((165542.487984203*rdxUxSq[0]+50077.05294843137*rdxLx[0]*rdxUx[0]+165542.487984203*rdxLxSq[0])*rho[0]-285331.1218356665*rho[0]*rdxLySq[1])*rdxUy[1]-4323.198815691917*rho[0]*rdxLyCu[1]+((-169865.6867998949*rdxUx[0])-169865.6867998949*rdxLx[0])*rho[0]*rdxLySq[1]+((-165542.487984203*rdxUxSq[0])-50077.05294843137*rdxLx[0]*rdxUx[0]-165542.487984203*rdxLxSq[0])*rho[0]*rdxLy[1])*volFac+((47596.75619199274*rdxUx[0]-47596.75619199274*rdxLx[0])*rdxUyCu[1]+((4001.037365484106*rdxUx[0]-4001.037365484106*rdxLx[0])*rdxLy[1]+46869.29485281381*rdxUxSq[0]-46869.29485281381*rdxLxSq[0])*rdxUySq[1]+((90880.70587313897*rdxLx[0]-90880.70587313897*rdxUx[0])*rdxLySq[1]+(91608.1672123179*rdxLxSq[0]-91608.1672123179*rdxUxSq[0])*rdxLy[1]-727.4613391789284*rdxUxCu[0]-48012.44838580926*rdxLx[0]*rdxUxSq[0]+48012.44838580926*rdxLxSq[0]*rdxUx[0]+727.4613391789284*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(1870.614872174387*rdxUx[0]*rdxUyCu[1]+(127201.8113078583*rdxUx[0]*rdxLy[1]+53728.21605078656*rdxUxSq[0]+97012.16573193281*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(127201.8113078583*rdxUx[0]*rdxLySq[1]+(276331.3858395386*rdxUxSq[0]+268329.3111085704*rdxLx[0]*rdxUx[0])*rdxLy[1]+53520.3699538783*rdxUxCu[0]+307144.569706189*rdxLx[0]*rdxUxSq[0]+98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+1870.614872174387*rdxUx[0]*rdxLyCu[1]+(53728.21605078656*rdxUxSq[0]+97012.16573193281*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(53520.3699538783*rdxUxCu[0]+307144.569706189*rdxLx[0]*rdxUxSq[0]+98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+1662.768775266122*rdxUxR4[0]+114523.1993964541*rdxLx[0]*rdxUxCu[0]+210548.0961680727*rdxLxSq[0]*rdxUxSq[0]+3117.691453623978*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((90880.70587313897*rdxLx[0]-90880.70587313897*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((4001.037365484106*rdxUx[0]-4001.037365484106*rdxLx[0])*rdxLySq[1]+(91608.1672123179*rdxLxSq[0]-91608.1672123179*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(47596.75619199274*rdxUx[0]-47596.75619199274*rdxLx[0])*rdxLyCu[1]+(46869.29485281381*rdxUxSq[0]-46869.29485281381*rdxLxSq[0])*rdxLySq[1]+((-727.4613391789284*rdxUxCu[0])-48012.44838580926*rdxLx[0]*rdxUxSq[0]+48012.44838580926*rdxLxSq[0]*rdxUx[0]+727.4613391789284*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-1870.614872174387*rdxLx[0]*rdxUyCu[1])+((-127201.8113078583*rdxLx[0]*rdxLy[1])-97012.16573193281*rdxLx[0]*rdxUx[0]-53728.21605078656*rdxLxSq[0])*rdxUySq[1]+((-127201.8113078583*rdxLx[0]*rdxLySq[1])+((-268329.3111085704*rdxLx[0]*rdxUx[0])-276331.3858395386*rdxLxSq[0])*rdxLy[1]-98259.24231338239*rdxLx[0]*rdxUxSq[0]-307144.569706189*rdxLxSq[0]*rdxUx[0]-53520.3699538783*rdxLxCu[0])*rdxUy[1]-1870.614872174387*rdxLx[0]*rdxLyCu[1]+((-97012.16573193281*rdxLx[0]*rdxUx[0])-53728.21605078656*rdxLxSq[0])*rdxLySq[1]+((-98259.24231338239*rdxLx[0]*rdxUxSq[0])-307144.569706189*rdxLxSq[0]*rdxUx[0]-53520.3699538783*rdxLxCu[0])*rdxLy[1]-3117.691453623978*rdxLx[0]*rdxUxCu[0]-210548.0961680727*rdxLxSq[0]*rdxUxSq[0]-114523.1993964541*rdxLxCu[0]*rdxUx[0]-1662.768775266122*rdxLxR4[0])*phiLx[3]+((-1584.0*rdxUyR4[1])+((-103032.0*rdxLy[1])-61482.0*rdxUx[0]-61482.0*rdxLx[0])*rdxUyCu[1]+(205848.0*rdxLySq[1]+(110061.0*rdxUx[0]+110061.0*rdxLx[0])*rdxLy[1]-58932.0*rdxUxSq[0]+97026.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*rdxUySq[1]+(3096.0*rdxLyCu[1]+(122403.0*rdxUx[0]+122403.0*rdxLx[0])*rdxLySq[1]+(120273.0*rdxUxSq[0]+151236.0*rdxLx[0]*rdxUx[0]+120273.0*rdxLxSq[0])*rdxLy[1]+966.0*rdxUxCu[0]+65688.0*rdxLx[0]*rdxUxSq[0]+65688.0*rdxLxSq[0]*rdxUx[0]+966.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-1944.0*rdxUx[0]*rdxUyCu[1])+((-132192.0*rdxUx[0]*rdxLy[1])-51552.0*rdxUxSq[0]-105102.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-132192.0*rdxUx[0]*rdxLySq[1])+((-287964.0*rdxUxSq[0])-278064.0*rdxLx[0]*rdxUx[0])*rdxLy[1]-51192.0*rdxUxCu[0]-319194.0*rdxLx[0]*rdxUxSq[0]-106542.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-1944.0*rdxUx[0]*rdxLyCu[1]+((-51552.0*rdxUxSq[0])-105102.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-51192.0*rdxUxCu[0])-319194.0*rdxLx[0]*rdxUxSq[0]-106542.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-1584.0*rdxUxR4[0]-109512.0*rdxLx[0]*rdxUxCu[0]-228312.0*rdxLxSq[0]*rdxUxSq[0]-3384.0*rdxLxCu[0]*rdxUx[0])*phiUx[2]+(3096.0*rdxLy[1]*rdxUyCu[1]+(205848.0*rdxLySq[1]+(122403.0*rdxUx[0]+122403.0*rdxLx[0])*rdxLy[1])*rdxUySq[1]+((-103032.0*rdxLyCu[1])+(110061.0*rdxUx[0]+110061.0*rdxLx[0])*rdxLySq[1]+(120273.0*rdxUxSq[0]+151236.0*rdxLx[0]*rdxUx[0]+120273.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-1584.0*rdxLyR4[1]+((-61482.0*rdxUx[0])-61482.0*rdxLx[0])*rdxLyCu[1]+((-58932.0*rdxUxSq[0])+97026.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*rdxLySq[1]+(966.0*rdxUxCu[0]+65688.0*rdxLx[0]*rdxUxSq[0]+65688.0*rdxLxSq[0]*rdxUx[0]+966.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-1944.0*rdxLx[0]*rdxUyCu[1])+((-132192.0*rdxLx[0]*rdxLy[1])-105102.0*rdxLx[0]*rdxUx[0]-51552.0*rdxLxSq[0])*rdxUySq[1]+((-132192.0*rdxLx[0]*rdxLySq[1])+((-278064.0*rdxLx[0]*rdxUx[0])-287964.0*rdxLxSq[0])*rdxLy[1]-106542.0*rdxLx[0]*rdxUxSq[0]-319194.0*rdxLxSq[0]*rdxUx[0]-51192.0*rdxLxCu[0])*rdxUy[1]-1944.0*rdxLx[0]*rdxLyCu[1]+((-105102.0*rdxLx[0]*rdxUx[0])-51552.0*rdxLxSq[0])*rdxLySq[1]+((-106542.0*rdxLx[0]*rdxUxSq[0])-319194.0*rdxLxSq[0]*rdxUx[0]-51192.0*rdxLxCu[0])*rdxLy[1]-3384.0*rdxLx[0]*rdxUxCu[0]-228312.0*rdxLxSq[0]*rdxUxSq[0]-109512.0*rdxLxCu[0]*rdxUx[0]-1584.0*rdxLxR4[0])*phiLx[2]+1496.491897739509*phiUy[0]*rdxUyR4[1]+((96897.85037863323*phiUy[0]+3367.106769913895*phiLy[0])*rdxLy[1]+(45180.0*rdxLx[0]-45180.0*rdxUx[0])*phiUy[1]-2340.0*rdxUx[0]*phiUx[1]+2340.0*rdxLx[0]*phiLx[1]+(57864.35337926103*phiUy[0]+2431.799333826703*phiUx[0])*rdxUx[0]+(57864.35337926103*phiUy[0]+2431.799333826703*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((224099.6616864915*phiLy[0]-224099.6616864915*phiUy[0])*rdxLySq[1]+((4950.0*rdxLx[0]-4950.0*rdxUx[0])*phiUy[1]-154440.0*rdxUx[0]*phiUx[1]+(98730.0*rdxLx[0]-98730.0*rdxUx[0])*phiLy[1]+154440.0*rdxLx[0]*phiLx[1]+((-136165.1742370272*phiUy[0])+160498.7560325623*phiUx[0]+133234.5442706207*phiLy[0])*rdxUx[0]+((-136165.1742370272*phiUy[0])+133234.5442706207*phiLy[0]+160498.7560325623*phiLx[0])*rdxLx[0])*rdxLy[1]+(44280.0*rdxLxSq[0]-44280.0*rdxUxSq[0])*phiUy[1]+((-63960.0*rdxUxSq[0])-119925.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(119925.0*rdxLx[0]*rdxUx[0]+63960.0*rdxLxSq[0])*phiLx[1]+(55172.74642429901*phiUy[0]+60930.08330865796*phiUx[0])*rdxUxSq[0]+((-125403.9425696018*phiUy[0])+130168.8143412238*phiUx[0]+130168.8143412238*phiLx[0])*rdxLx[0]*rdxUx[0]+(55172.74642429901*phiUy[0]+60930.08330865796*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-3367.106769913895*phiUy[0])-96897.85037863323*phiLy[0])*rdxLyCu[1]+((98730.0*rdxUx[0]-98730.0*rdxLx[0])*phiUy[1]+154440.0*rdxUx[0]*phiUx[1]+(4950.0*rdxUx[0]-4950.0*rdxLx[0])*phiLy[1]-154440.0*rdxLx[0]*phiLx[1]+((-133234.5442706207*phiUy[0])-160498.7560325623*phiUx[0]+136165.1742370272*phiLy[0])*rdxUx[0]+((-133234.5442706207*phiUy[0])+136165.1742370272*phiLy[0]-160498.7560325623*phiLx[0])*rdxLx[0])*rdxLySq[1]+((99630.0*rdxUxSq[0]-99630.0*rdxLxSq[0])*phiUy[1]+(99630.0*rdxLxSq[0]-99630.0*rdxUxSq[0])*phiLy[1]+(131062.5525579294*phiLy[0]-131062.5525579294*phiUy[0])*rdxUxSq[0]+(181740.6271365871*phiLy[0]-181740.6271365871*phiUy[0])*rdxLx[0]*rdxUx[0]+(131062.5525579294*phiLy[0]-131062.5525579294*phiUy[0])*rdxLxSq[0])*rdxLy[1]+(900.0*rdxUxCu[0]+59400.0*rdxLx[0]*rdxUxSq[0]-59400.0*rdxLxSq[0]*rdxUx[0]-900.0*rdxLxCu[0])*phiUy[1]+((-61620.0*rdxUxCu[0])-27105.0*rdxLx[0]*rdxUxSq[0]-117585.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(117585.0*rdxLx[0]*rdxUxSq[0]+27105.0*rdxLxSq[0]*rdxUx[0]+61620.0*rdxLxCu[0])*phiLx[1]+(58498.28397483125*phiUx[0]-1195.115057222525*phiUy[0])*rdxUxCu[0]+((-81267.82389113172*phiUy[0])+28168.34228349264*phiUx[0]+127737.0150073971*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-81267.82389113172*phiUy[0])+127737.0150073971*phiUx[0]+28168.34228349264*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(58498.28397483125*phiLx[0]-1195.115057222525*phiUy[0])*rdxLxCu[0])*rdxUy[1]-1496.491897739509*phiLy[0]*rdxLyR4[1]+(2340.0*rdxUx[0]*phiUx[1]+(45180.0*rdxUx[0]-45180.0*rdxLx[0])*phiLy[1]-2340.0*rdxLx[0]*phiLx[1]+((-2431.799333826703*phiUx[0])-57864.35337926103*phiLy[0])*rdxUx[0]+((-57864.35337926103*phiLy[0])-2431.799333826703*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((63960.0*rdxUxSq[0]+119925.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(44280.0*rdxUxSq[0]-44280.0*rdxLxSq[0])*phiLy[1]+((-119925.0*rdxLx[0]*rdxUx[0])-63960.0*rdxLxSq[0])*phiLx[1]+((-60930.08330865796*phiUx[0])-55172.74642429901*phiLy[0])*rdxUxSq[0]+((-130168.8143412238*phiUx[0])+125403.9425696018*phiLy[0]-130168.8143412238*phiLx[0])*rdxLx[0]*rdxUx[0]+((-55172.74642429901*phiLy[0])-60930.08330865796*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((61620.0*rdxUxCu[0]+27105.0*rdxLx[0]*rdxUxSq[0]+117585.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-900.0*rdxUxCu[0])-59400.0*rdxLx[0]*rdxUxSq[0]+59400.0*rdxLxSq[0]*rdxUx[0]+900.0*rdxLxCu[0])*phiLy[1]+((-117585.0*rdxLx[0]*rdxUxSq[0])-27105.0*rdxLxSq[0]*rdxUx[0]-61620.0*rdxLxCu[0])*phiLx[1]+(1195.115057222525*phiLy[0]-58498.28397483125*phiUx[0])*rdxUxCu[0]+((-28168.34228349264*phiUx[0])+81267.82389113172*phiLy[0]-127737.0150073971*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-127737.0150073971*phiUx[0])+81267.82389113172*phiLy[0]-28168.34228349264*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(1195.115057222525*phiLy[0]-58498.28397483125*phiLx[0])*rdxLxCu[0])*rdxLy[1]))/(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0]); 
  phiC[3] = (((288.0*rdxUyCu[1]+(19584.0*rdxLy[1]+9224.0*rdxUx[0]+9224.0*rdxLx[0])*rdxUySq[1]+(19584.0*rdxLySq[1]+(42368.0*rdxUx[0]+42368.0*rdxLx[0])*rdxLy[1]+9224.0*rdxUxSq[0]+42368.0*rdxLx[0]*rdxUx[0]+9224.0*rdxLxSq[0])*rdxUy[1]+288.0*rdxLyCu[1]+(9224.0*rdxUx[0]+9224.0*rdxLx[0])*rdxLySq[1]+(9224.0*rdxUxSq[0]+42368.0*rdxLx[0]*rdxUx[0]+9224.0*rdxLxSq[0])*rdxLy[1]+288.0*rdxUxCu[0]+19584.0*rdxLx[0]*rdxUxSq[0]+19584.0*rdxLxSq[0]*rdxUx[0]+288.0*rdxLxCu[0])*rho[3]+((21435.86079447242*rdxLx[0]-21435.86079447242*rdxUx[0])*rdxUySq[1]+((3962.932247717591*rdxUx[0]-3962.932247717591*rdxLx[0])*rdxLy[1]-22156.39393042108*rdxUxSq[0]+22156.39393042108*rdxLxSq[0])*rdxUy[1]+(21435.86079447242*rdxLx[0]-21435.86079447242*rdxUx[0])*rdxLySq[1]+(22156.39393042108*rdxLxSq[0]-22156.39393042108*rdxUxSq[0])*rdxLy[1]-720.5331359486529*rdxUxCu[0]-47555.18697261109*rdxLx[0]*rdxUxSq[0]+47555.18697261109*rdxLxSq[0]*rdxUx[0]+720.5331359486529*rdxLxCu[0])*rho[2]+((-720.5331359486529*rdxUyCu[1])+((-47555.18697261109*rdxLy[1])-22156.39393042108*rdxUx[0]-22156.39393042108*rdxLx[0])*rdxUySq[1]+(47555.18697261109*rdxLySq[1]-21435.86079447242*rdxUxSq[0]+3962.932247717591*rdxLx[0]*rdxUx[0]-21435.86079447242*rdxLxSq[0])*rdxUy[1]+720.5331359486529*rdxLyCu[1]+(22156.39393042108*rdxUx[0]+22156.39393042108*rdxLx[0])*rdxLySq[1]+(21435.86079447242*rdxUxSq[0]-3962.932247717591*rdxLx[0]*rdxUx[0]+21435.86079447242*rdxLxSq[0])*rdxLy[1])*rho[1]+(55432.0*rdxUx[0]-55432.0*rdxLx[0])*rho[0]*rdxUySq[1]+(55432.0*rdxUxSq[0]-55432.0*rdxLxSq[0])*rho[0]*rdxUy[1]+(55432.0*rdxLx[0]-55432.0*rdxUx[0])*rho[0]*rdxLySq[1]+(55432.0*rdxLxSq[0]-55432.0*rdxUxSq[0])*rho[0]*rdxLy[1])*volFac+(528.0*rdxUyR4[1]+(34344.0*rdxLy[1]+15914.0*rdxUx[0]+15914.0*rdxLx[0])*rdxUyCu[1]+((-68616.0*rdxLySq[1])+((-37072.0*rdxUx[0])-37072.0*rdxLx[0])*rdxLy[1]+15134.0*rdxUxSq[0]-41362.0*rdxLx[0]*rdxUx[0]+15134.0*rdxLxSq[0])*rdxUySq[1]+((-1032.0*rdxLyCu[1])+((-32056.0*rdxUx[0])-32056.0*rdxLx[0])*rdxLySq[1]+((-31276.0*rdxUxSq[0])-32782.0*rdxLx[0]*rdxUx[0]-31276.0*rdxLxSq[0])*rdxLy[1]-252.0*rdxUxCu[0]-17136.0*rdxLx[0]*rdxUxSq[0]-17136.0*rdxLxSq[0]*rdxUx[0]-252.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-252.0*rdxUx[0]*rdxUyCu[1])+((-17136.0*rdxUx[0]*rdxLy[1])+15134.0*rdxUxSq[0]-31276.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-17136.0*rdxUx[0]*rdxLySq[1])+((-41362.0*rdxUxSq[0])-32782.0*rdxLx[0]*rdxUx[0])*rdxLy[1]+15914.0*rdxUxCu[0]-37072.0*rdxLx[0]*rdxUxSq[0]-32056.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-252.0*rdxUx[0]*rdxLyCu[1]+(15134.0*rdxUxSq[0]-31276.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(15914.0*rdxUxCu[0]-37072.0*rdxLx[0]*rdxUxSq[0]-32056.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+528.0*rdxUxR4[0]+34344.0*rdxLx[0]*rdxUxCu[0]-68616.0*rdxLxSq[0]*rdxUxSq[0]-1032.0*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((-1032.0*rdxLy[1]*rdxUyCu[1])+(((-32056.0*rdxUx[0])-32056.0*rdxLx[0])*rdxLy[1]-68616.0*rdxLySq[1])*rdxUySq[1]+(34344.0*rdxLyCu[1]+((-37072.0*rdxUx[0])-37072.0*rdxLx[0])*rdxLySq[1]+((-31276.0*rdxUxSq[0])-32782.0*rdxLx[0]*rdxUx[0]-31276.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+528.0*rdxLyR4[1]+(15914.0*rdxUx[0]+15914.0*rdxLx[0])*rdxLyCu[1]+(15134.0*rdxUxSq[0]-41362.0*rdxLx[0]*rdxUx[0]+15134.0*rdxLxSq[0])*rdxLySq[1]+((-252.0*rdxUxCu[0])-17136.0*rdxLx[0]*rdxUxSq[0]-17136.0*rdxLxSq[0]*rdxUx[0]-252.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-252.0*rdxLx[0]*rdxUyCu[1])+((-17136.0*rdxLx[0]*rdxLy[1])-31276.0*rdxLx[0]*rdxUx[0]+15134.0*rdxLxSq[0])*rdxUySq[1]+((-17136.0*rdxLx[0]*rdxLySq[1])+((-32782.0*rdxLx[0]*rdxUx[0])-41362.0*rdxLxSq[0])*rdxLy[1]-32056.0*rdxLx[0]*rdxUxSq[0]-37072.0*rdxLxSq[0]*rdxUx[0]+15914.0*rdxLxCu[0])*rdxUy[1]-252.0*rdxLx[0]*rdxLyCu[1]+(15134.0*rdxLxSq[0]-31276.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-32056.0*rdxLx[0]*rdxUxSq[0])-37072.0*rdxLxSq[0]*rdxUx[0]+15914.0*rdxLxCu[0])*rdxLy[1]-1032.0*rdxLx[0]*rdxUxCu[0]-68616.0*rdxLxSq[0]*rdxUxSq[0]+34344.0*rdxLxCu[0]*rdxUx[0]+528.0*rdxLxR4[0])*phiLx[3]+((20625.26101653019*rdxLx[0]-20625.26101653019*rdxUx[0])*rdxUyCu[1]+((1733.782858376446*rdxLx[0]-1733.782858376446*rdxUx[0])*rdxLy[1]-20310.02776955265*rdxUxSq[0]+20310.02776955265*rdxLxSq[0])*rdxUySq[1]+((39381.63921169356*rdxUx[0]-39381.63921169356*rdxLx[0])*rdxLySq[1]+(39696.8724586711*rdxUxSq[0]-39696.8724586711*rdxLxSq[0])*rdxLy[1]+315.2332469775357*rdxUxCu[0]+20805.39430051735*rdxLx[0]*rdxUxSq[0]-20805.39430051735*rdxLxSq[0]*rdxUx[0]-315.2332469775357*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(311.7691453623978*rdxUx[0]*rdxUyCu[1]+(21200.30188464305*rdxUx[0]*rdxLy[1]-14130.0704881469*rdxUxSq[0]+34100.61629941605*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(21200.30188464305*rdxUx[0]*rdxLySq[1]+(50323.00416310616*rdxUxSq[0]+41406.40660574158*rdxLx[0]*rdxUx[0])*rdxLy[1]-14940.67026608914*rdxUxCu[0]+45864.70538442387*rdxLx[0]*rdxUxSq[0]+34911.21607735829*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+311.7691453623978*rdxUx[0]*rdxLyCu[1]+(34100.61629941605*rdxLx[0]*rdxUx[0]-14130.0704881469*rdxUxSq[0])*rdxLySq[1]+((-14940.67026608914*rdxUxCu[0])+45864.70538442387*rdxLx[0]*rdxUxSq[0]+34911.21607735829*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-498.8306325798365*rdxUxR4[0]-32299.28345954441*rdxLx[0]*rdxUxCu[0]+74699.88722883051*rdxLxSq[0]*rdxUxSq[0]+1122.368923304632*rdxLxCu[0]*rdxUx[0])*phiUx[2]+((39381.63921169356*rdxUx[0]-39381.63921169356*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((1733.782858376446*rdxLx[0]-1733.782858376446*rdxUx[0])*rdxLySq[1]+(39696.8724586711*rdxUxSq[0]-39696.8724586711*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(20625.26101653019*rdxLx[0]-20625.26101653019*rdxUx[0])*rdxLyCu[1]+(20310.02776955265*rdxLxSq[0]-20310.02776955265*rdxUxSq[0])*rdxLySq[1]+(315.2332469775357*rdxUxCu[0]+20805.39430051735*rdxLx[0]*rdxUxSq[0]-20805.39430051735*rdxLxSq[0]*rdxUx[0]-315.2332469775357*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-311.7691453623978*rdxLx[0]*rdxUyCu[1])+((-21200.30188464305*rdxLx[0]*rdxLy[1])-34100.61629941605*rdxLx[0]*rdxUx[0]+14130.0704881469*rdxLxSq[0])*rdxUySq[1]+((-21200.30188464305*rdxLx[0]*rdxLySq[1])+((-41406.40660574158*rdxLx[0]*rdxUx[0])-50323.00416310616*rdxLxSq[0])*rdxLy[1]-34911.21607735829*rdxLx[0]*rdxUxSq[0]-45864.70538442387*rdxLxSq[0]*rdxUx[0]+14940.67026608914*rdxLxCu[0])*rdxUy[1]-311.7691453623978*rdxLx[0]*rdxLyCu[1]+(14130.0704881469*rdxLxSq[0]-34100.61629941605*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-34911.21607735829*rdxLx[0]*rdxUxSq[0])-45864.70538442387*rdxLxSq[0]*rdxUx[0]+14940.67026608914*rdxLxCu[0])*rdxLy[1]-1122.368923304632*rdxLx[0]*rdxUxCu[0]-74699.88722883051*rdxLxSq[0]*rdxUxSq[0]+32299.28345954441*rdxLxCu[0]*rdxUx[0]+498.8306325798365*rdxLxR4[0])*phiLx[2]-498.8306325798365*phiUy[1]*rdxUyR4[1]+(((-32299.28345954441*phiUy[1])-1122.368923304632*phiLy[1])*rdxLy[1]+((-14940.67026608914*rdxUx[0])-14940.67026608914*rdxLx[0])*phiUy[1]+315.2332469775357*rdxUx[0]*phiUx[1]+315.2332469775357*rdxLx[0]*phiLx[1]+(19578.0*phiUy[0]-390.0*phiUx[0])*rdxUx[0]+(390.0*phiLx[0]-19578.0*phiUy[0])*rdxLx[0])*rdxUyCu[1]+((74699.88722883051*phiUy[1]-74699.88722883051*phiLy[1])*rdxLySq[1]+((45864.70538442387*rdxUx[0]+45864.70538442387*rdxLx[0])*phiUy[1]+20805.39430051735*rdxUx[0]*phiUx[1]+((-34911.21607735829*rdxUx[0])-34911.21607735829*rdxLx[0])*phiLy[1]+20805.39430051735*rdxLx[0]*phiLx[1]+(2145.0*phiUy[0]-25740.0*phiUx[0]+42783.0*phiLy[0])*rdxUx[0]+((-2145.0*phiUy[0])-42783.0*phiLy[0]+25740.0*phiLx[0])*rdxLx[0])*rdxLy[1]+((-14130.0704881469*rdxUxSq[0])+50323.00416310616*rdxLx[0]*rdxUx[0]-14130.0704881469*rdxLxSq[0])*phiUy[1]+(39696.8724586711*rdxLx[0]*rdxUx[0]-20310.02776955265*rdxUxSq[0])*phiUx[1]+(39696.8724586711*rdxLx[0]*rdxUx[0]-20310.02776955265*rdxLxSq[0])*phiLx[1]+(19188.0*phiUy[0]+19188.0*phiUx[0])*rdxUxSq[0]+(43173.0*phiLx[0]-43173.0*phiUx[0])*rdxLx[0]*rdxUx[0]+((-19188.0*phiUy[0])-19188.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((1122.368923304632*phiUy[1]+32299.28345954441*phiLy[1])*rdxLyCu[1]+((34911.21607735829*rdxUx[0]+34911.21607735829*rdxLx[0])*phiUy[1]-20805.39430051735*rdxUx[0]*phiUx[1]+((-45864.70538442387*rdxUx[0])-45864.70538442387*rdxLx[0])*phiLy[1]-20805.39430051735*rdxLx[0]*phiLx[1]+((-42783.0*phiUy[0])+25740.0*phiUx[0]-2145.0*phiLy[0])*rdxUx[0]+(42783.0*phiUy[0]+2145.0*phiLy[0]-25740.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((34100.61629941605*rdxUxSq[0]+41406.40660574158*rdxLx[0]*rdxUx[0]+34100.61629941605*rdxLxSq[0])*phiUy[1]+((-34100.61629941605*rdxUxSq[0])-41406.40660574158*rdxLx[0]*rdxUx[0]-34100.61629941605*rdxLxSq[0])*phiLy[1]+(43173.0*phiLy[0]-43173.0*phiUy[0])*rdxUxSq[0]+(43173.0*phiUy[0]-43173.0*phiLy[0])*rdxLxSq[0])*rdxLy[1]+(311.7691453623978*rdxUxCu[0]+21200.30188464305*rdxLx[0]*rdxUxSq[0]+21200.30188464305*rdxLxSq[0]*rdxUx[0]+311.7691453623978*rdxLxCu[0])*phiUy[1]+((-20625.26101653019*rdxUxCu[0])-1733.782858376446*rdxLx[0]*rdxUxSq[0]+39381.63921169356*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(39381.63921169356*rdxLx[0]*rdxUxSq[0]-1733.782858376446*rdxLxSq[0]*rdxUx[0]-20625.26101653019*rdxLxCu[0])*phiLx[1]+(19578.0*phiUx[0]-390.0*phiUy[0])*rdxUxCu[0]+((-25740.0*phiUy[0])+2145.0*phiUx[0]+42783.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(25740.0*phiUy[0]-42783.0*phiUx[0]-2145.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(390.0*phiUy[0]-19578.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+498.8306325798365*phiLy[1]*rdxLyR4[1]+((-315.2332469775357*rdxUx[0]*phiUx[1])+(14940.67026608914*rdxUx[0]+14940.67026608914*rdxLx[0])*phiLy[1]-315.2332469775357*rdxLx[0]*phiLx[1]+(390.0*phiUx[0]-19578.0*phiLy[0])*rdxUx[0]+(19578.0*phiLy[0]-390.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((20310.02776955265*rdxUxSq[0]-39696.8724586711*rdxLx[0]*rdxUx[0])*phiUx[1]+(14130.0704881469*rdxUxSq[0]-50323.00416310616*rdxLx[0]*rdxUx[0]+14130.0704881469*rdxLxSq[0])*phiLy[1]+(20310.02776955265*rdxLxSq[0]-39696.8724586711*rdxLx[0]*rdxUx[0])*phiLx[1]+((-19188.0*phiUx[0])-19188.0*phiLy[0])*rdxUxSq[0]+(43173.0*phiUx[0]-43173.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(19188.0*phiLy[0]+19188.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((20625.26101653019*rdxUxCu[0]+1733.782858376446*rdxLx[0]*rdxUxSq[0]-39381.63921169356*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-311.7691453623978*rdxUxCu[0])-21200.30188464305*rdxLx[0]*rdxUxSq[0]-21200.30188464305*rdxLxSq[0]*rdxUx[0]-311.7691453623978*rdxLxCu[0])*phiLy[1]+((-39381.63921169356*rdxLx[0]*rdxUxSq[0])+1733.782858376446*rdxLxSq[0]*rdxUx[0]+20625.26101653019*rdxLxCu[0])*phiLx[1]+(390.0*phiLy[0]-19578.0*phiUx[0])*rdxUxCu[0]+((-2145.0*phiUx[0])+25740.0*phiLy[0]-42783.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(42783.0*phiUx[0]-25740.0*phiLy[0]+2145.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(19578.0*phiLx[0]-390.0*phiLy[0])*rdxLxCu[0])*rdxLy[1])/(48.0*rdxUyR4[1]+(6432.0*rdxLy[1]+3362.0*rdxUx[0]+3362.0*rdxLx[0])*rdxUyCu[1]+(215568.0*rdxLySq[1]+(228616.0*rdxUx[0]+228616.0*rdxLx[0])*rdxLy[1]+6628.0*rdxUxSq[0]+225546.0*rdxLx[0]*rdxUx[0]+6628.0*rdxLxSq[0])*rdxUySq[1]+(6432.0*rdxLyCu[1]+(228616.0*rdxUx[0]+228616.0*rdxLx[0])*rdxLySq[1]+(225546.0*rdxUxSq[0]+470072.0*rdxLx[0]*rdxUx[0]+225546.0*rdxLxSq[0])*rdxLy[1]+3362.0*rdxUxCu[0]+228616.0*rdxLx[0]*rdxUxSq[0]+228616.0*rdxLxSq[0]*rdxUx[0]+3362.0*rdxLxCu[0])*rdxUy[1]+48.0*rdxLyR4[1]+(3362.0*rdxUx[0]+3362.0*rdxLx[0])*rdxLyCu[1]+(6628.0*rdxUxSq[0]+225546.0*rdxLx[0]*rdxUx[0]+6628.0*rdxLxSq[0])*rdxLySq[1]+(3362.0*rdxUxCu[0]+228616.0*rdxLx[0]*rdxUxSq[0]+228616.0*rdxLxSq[0]*rdxUx[0]+3362.0*rdxLxCu[0])*rdxLy[1]+48.0*rdxUxR4[0]+6432.0*rdxLx[0]*rdxUxCu[0]+215568.0*rdxLxSq[0]*rdxUxSq[0]+6432.0*rdxLxCu[0]*rdxUx[0]+48.0*rdxLxR4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((1.732050807568877*rdxCp2[0]*rho[1]+6.0*rho[0]*rdxCp2[1]+50.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+5.196152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]+9.0*phiUx[0]+75.0*phiLy[0]+72.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiUx[1]+(90.0*phiUx[0]+480.0*bcVals[0])*rdxCp2Sq[0])/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[1]+17.32050807568877*rdxCp2[0]*rho[0])*volFac+((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+27.0*rdxCp2[0]*phiLy[1]+(25.98076211353316*phiUy[0]+51.96152422706631*phiUx[0]+25.98076211353316*phiLy[0]-415.6921938165305*bcVals[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiUx[1]+(103.9230484541326*phiUx[0]-415.6921938165305*bcVals[0])*rdxCp2Sq[0])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = ((1.732050807568877*rdxCp2[0]*rho[3]+(40.0*rdxCp2[1]+50.0*rdxCp2[0])*rho[2])*volFac-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUy[3]+((-138.5640646055102*rdxCp2[0]*rdxCp2[1])-207.8460969082653*rdxCp2Sq[0])*phiUx[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-400.0*rdxCp2Sq[1])-500.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(120.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*phiUx[2]+((-400.0*rdxCp2Sq[1])-500.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(346.4101615137754*phiUy[0]-346.4101615137754*phiLy[0])*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUy[1]-30.0*rdxCp2[0]*phiLy[1]+(433.0127018922193*phiUy[0]-433.0127018922193*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[3] = (((40.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[3]+17.32050807568877*rdxCp2[0]*rho[2])*volFac+((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-800.0*rdxCp2[0]*rdxCp2[1])-240.0*rdxCp2Sq[0])*phiUx[3]+((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(692.8203230275509*rdxCp2[0]*rdxCp2[1]+207.8460969082653*rdxCp2Sq[0])*phiUx[2]-173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(692.8203230275509*phiUy[1]-692.8203230275509*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(150.0*phiUy[0]-150.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(3200.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+840.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((10.39230484541326*rdxCp2[0]*rho[1]-18.0*rho[0]*rdxCp2[1]-60.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-27.0*phiUy[0])-27.0*phiLy[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+41.56921938165305*rdxCp2[0]*phiUx[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+((-90.0*phiUy[0])-36.0*phiUx[0]-90.0*phiLy[0]+24.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiUx[1]+(160.0*bcVals[0]-60.0*phiUx[0])*rdxCp2Sq[0]))/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[1]-17.32050807568877*rdxCp2[0]*rho[0])*volFac+((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+18.0*rdxCp2[0]*phiLy[1]+((-25.98076211353316*phiUy[0])+51.96152422706631*phiUx[0]-25.98076211353316*phiLy[0]+69.28203230275508*bcVals[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*bcVals[0]*rdxCp2Sq[0])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((5.196152422706631*rdxCp2[0]*rho[3]+((-60.0*rdxCp2[1])-30.0*rdxCp2[0])*rho[2])*volFac-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUy[3]+(277.1281292110203*rdxCp2[0]*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0])*phiUx[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(600.0*rdxCp2Sq[1]+300.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-240.0*rdxCp2[0]*rdxCp2[1])-60.0*rdxCp2Sq[0])*phiUx[2]+(600.0*rdxCp2Sq[1]+300.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(519.6152422706631*phiLy[0]-519.6152422706631*phiUy[0])*rdxCp2Sq[1]+(90.0*rdxCp2[0]*phiUy[1]-90.0*rdxCp2[0]*phiLy[1]+(259.8076211353315*phiLy[0]-259.8076211353315*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[3] = (((30.0*rdxCp2[1]+3.0*rdxCp2[0])*rho[3]-8.660254037844386*rdxCp2[0]*rho[2])*volFac+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[2]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiLy[0]-75.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((1.732050807568877*rdxCp2[0]*rho[1]-6.0*rho[0]*rdxCp2[1]-50.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-9.0*phiUy[0])-9.0*phiLy[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]+5.196152422706631*rdxCp2[0]*phiLy[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-72.0*rdxCp2[0]*bcVals[1]+((-75.0*phiUy[0])-75.0*phiLy[0]-9.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiLx[1]-480.0*rdxCp2Sq[0]*bcVals[1]-90.0*phiLx[0]*rdxCp2Sq[0]))/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[1]-17.32050807568877*rdxCp2[0]*rho[0])*volFac+((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]+27.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]+415.6921938165305*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiUy[0])-25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiLx[1]+415.6921938165305*rdxCp2Sq[0]*bcVals[1]-103.9230484541326*phiLx[0]*rdxCp2Sq[0])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((1.732050807568877*rdxCp2[0]*rho[3]+((-40.0*rdxCp2[1])-50.0*rdxCp2[0])*rho[2])*volFac-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUy[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-138.5640646055102*rdxCp2[0]*rdxCp2[1])-207.8460969082653*rdxCp2Sq[0])*phiLx[3]+(400.0*rdxCp2Sq[1]+500.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(400.0*rdxCp2Sq[1]+500.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-120.0*rdxCp2[0]*rdxCp2[1])-180.0*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*phiLy[0]-346.4101615137754*phiUy[0])*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUy[1]-30.0*rdxCp2[0]*phiLy[1]+(433.0127018922193*phiLy[0]-433.0127018922193*phiUy[0])*rdxCp2[0])*rdxCp2[1]))/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[3] = (((40.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[3]-17.32050807568877*rdxCp2[0]*rho[2])*volFac+((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-800.0*rdxCp2[0]*rdxCp2[1])-240.0*rdxCp2Sq[0])*phiLx[3]+173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[2]+173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-692.8203230275509*rdxCp2[0]*rdxCp2[1])-207.8460969082653*rdxCp2Sq[0])*phiLx[2]+(692.8203230275509*phiUy[1]-692.8203230275509*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(150.0*phiLy[0]-150.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])/(3200.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+840.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((10.39230484541326*rdxCp2[0]*rho[1]+18.0*rho[0]*rdxCp2[1]+60.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*phiUy[0]+27.0*phiLy[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+41.56921938165305*rdxCp2[0]*phiLx[1]+24.0*rdxCp2[0]*bcVals[1]+(90.0*phiUy[0]+90.0*phiLy[0]+36.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiLx[1]+160.0*rdxCp2Sq[0]*bcVals[1]+60.0*phiLx[0]*rdxCp2Sq[0])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[1]+17.32050807568877*rdxCp2[0]*rho[0])*volFac+((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]+18.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiUy[0]+25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*bcVals[1])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = ((5.196152422706631*rdxCp2[0]*rho[3]+(60.0*rdxCp2[1]+30.0*rdxCp2[0])*rho[2])*volFac-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUy[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(277.1281292110203*rdxCp2[0]*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0])*phiLx[3]+((-600.0*rdxCp2Sq[1])-300.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-600.0*rdxCp2Sq[1])-300.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(240.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiLx[2]+(519.6152422706631*phiUy[0]-519.6152422706631*phiLy[0])*rdxCp2Sq[1]+(90.0*rdxCp2[0]*phiUy[1]-90.0*rdxCp2[0]*phiLy[1]+(259.8076211353315*phiUy[0]-259.8076211353315*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[3] = (((30.0*rdxCp2[1]+3.0*rdxCp2[0])*rho[3]+8.660254037844386*rdxCp2[0]*rho[2])*volFac+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]-346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]-75.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((1.732050807568877*rdxCp2[1]*rho[2]+50.0*rho[0]*rdxCp2[1]+6.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(480.0*rdxCp2Sq[1]+72.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+90.0*phiUy[0]*rdxCp2Sq[1]+((-86.60254037844386*rdxCp2[0]*phiUx[1])+86.60254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUy[0]+75.0*phiUx[0]+75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-10.39230484541326*rdxCp2Sq[0]*phiUx[1]+10.39230484541326*rdxCp2Sq[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0])*rdxCp2Sq[0])/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = ((1.732050807568877*rdxCp2[1]*rho[3]+(50.0*rdxCp2[1]+40.0*rdxCp2[0])*rho[1])*volFac+((-207.8460969082653*rdxCp2Sq[1])-138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiUy[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+180.0*phiUy[1]*rdxCp2Sq[1]+(120.0*rdxCp2[0]*phiUy[1]-500.0*rdxCp2[0]*phiUx[1]-500.0*rdxCp2[0]*phiLx[1]+(433.0127018922193*phiUx[0]-433.0127018922193*phiLx[0])*rdxCp2[0])*rdxCp2[1]-400.0*rdxCp2Sq[0]*phiUx[1]-400.0*rdxCp2Sq[0]*phiLx[1]+(346.4101615137754*phiUx[0]-346.4101615137754*phiLx[0])*rdxCp2Sq[0])/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 
  phiC[2] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[2]+17.32050807568877*rho[0]*rdxCp2[1])*volFac+((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]+((-415.6921938165305*rdxCp2Sq[1])-415.6921938165305*rdxCp2[0]*rdxCp2[1])*bcVals[2]+103.9230484541326*phiUy[0]*rdxCp2Sq[1]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]+25.98076211353316*phiUx[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((9.0*rdxCp2[1]+40.0*rdxCp2[0])*rho[3]+17.32050807568877*rdxCp2[1]*rho[1])*volFac+((-240.0*rdxCp2Sq[1])-800.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiUx[3]+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiLx[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-692.8203230275509*rdxCp2Sq[0])*phiLx[2]+207.8460969082653*phiUy[1]*rdxCp2Sq[1]+(692.8203230275509*rdxCp2[0]*phiUy[1]-173.2050807568877*rdxCp2[0]*phiUx[1]-173.2050807568877*rdxCp2[0]*phiLx[1]+(150.0*phiUx[0]-150.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(840.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+3200.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((10.39230484541326*rdxCp2[1]*rho[2]-60.0*rho[0]*rdxCp2[1]-18.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiUy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]-60.0*phiUy[0]*rdxCp2Sq[1]+(103.9230484541326*rdxCp2[0]*phiUx[1]-103.9230484541326*rdxCp2[0]*phiLx[1]+((-36.0*phiUy[0])-90.0*phiUx[0]-90.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0]*phiUx[1]-31.17691453623978*rdxCp2Sq[0]*phiLx[1]+((-27.0*phiUx[0])-27.0*phiLx[0])*rdxCp2Sq[0]))/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((5.196152422706631*rdxCp2[1]*rho[3]+((-30.0*rdxCp2[1])-60.0*rdxCp2[0])*rho[1])*volFac+(69.28203230275508*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*phiUy[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUx[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLx[3]+90.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-90.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-60.0*phiUy[1]*rdxCp2Sq[1]+((-240.0*rdxCp2[0]*phiUy[1])+300.0*rdxCp2[0]*phiUx[1]+300.0*rdxCp2[0]*phiLx[1]+(259.8076211353315*phiLx[0]-259.8076211353315*phiUx[0])*rdxCp2[0])*rdxCp2[1]+600.0*rdxCp2Sq[0]*phiUx[1]+600.0*rdxCp2Sq[0]*phiLx[1]+(519.6152422706631*phiLx[0]-519.6152422706631*phiUx[0])*rdxCp2Sq[0]))/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 
  phiC[2] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[2]-17.32050807568877*rho[0]*rdxCp2[1])*volFac+((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]-60.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]-25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((3.0*rdxCp2[1]+30.0*rdxCp2[0])*rho[3]-8.660254037844386*rdxCp2[1]*rho[1])*volFac-400.0*rdxCp2[0]*rdxCp2[1]*phiUy[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*rdxCp2[0]*phiUy[1]+86.60254037844386*rdxCp2[0]*phiUx[1]+86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiLx[0]-75.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((1.732050807568877*rdxCp2[1]*rho[2]-50.0*rho[0]*rdxCp2[1]-6.0*rdxCp2[0]*rho[0])*volFac-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-480.0*rdxCp2Sq[1])-72.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]-90.0*phiLy[0]*rdxCp2Sq[1]+(86.60254037844386*rdxCp2[0]*phiUx[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+((-75.0*phiUx[0])-9.0*phiLy[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]+10.39230484541326*rdxCp2Sq[0]*phiUx[1]-10.39230484541326*rdxCp2Sq[0]*phiLx[1]+((-9.0*phiUx[0])-9.0*phiLx[0])*rdxCp2Sq[0]))/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((1.732050807568877*rdxCp2[1]*rho[3]+((-50.0*rdxCp2[1])-40.0*rdxCp2[0])*rho[1])*volFac-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-207.8460969082653*rdxCp2Sq[1])-138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiLy[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]-180.0*phiLy[1]*rdxCp2Sq[1]+(500.0*rdxCp2[0]*phiUx[1]-120.0*rdxCp2[0]*phiLy[1]+500.0*rdxCp2[0]*phiLx[1]+(433.0127018922193*phiLx[0]-433.0127018922193*phiUx[0])*rdxCp2[0])*rdxCp2[1]+400.0*rdxCp2Sq[0]*phiUx[1]+400.0*rdxCp2Sq[0]*phiLx[1]+(346.4101615137754*phiLx[0]-346.4101615137754*phiUx[0])*rdxCp2Sq[0]))/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 
  phiC[2] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[2]-17.32050807568877*rho[0]*rdxCp2[1])*volFac+((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+(415.6921938165305*rdxCp2Sq[1]+415.6921938165305*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]-103.9230484541326*phiLy[0]*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+((-25.98076211353316*phiUx[0])-51.96152422706631*phiLy[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((9.0*rdxCp2[1]+40.0*rdxCp2[0])*rho[3]-17.32050807568877*rdxCp2[1]*rho[1])*volFac+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiUx[3]+((-240.0*rdxCp2Sq[1])-800.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiLx[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-692.8203230275509*rdxCp2Sq[0])*phiLx[2]-207.8460969082653*phiLy[1]*rdxCp2Sq[1]+(173.2050807568877*rdxCp2[0]*phiUx[1]-692.8203230275509*rdxCp2[0]*phiLy[1]+173.2050807568877*rdxCp2[0]*phiLx[1]+(150.0*phiLx[0]-150.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])/(840.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+3200.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((10.39230484541326*rdxCp2[1]*rho[2]+60.0*rho[0]*rdxCp2[1]+18.0*rdxCp2[0]*rho[0])*volFac-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiLy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+60.0*phiLy[0]*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2[0]*phiUx[1])+103.9230484541326*rdxCp2[0]*phiLx[1]+(90.0*phiUx[0]+36.0*phiLy[0]+90.0*phiLx[0])*rdxCp2[0])*rdxCp2[1]-31.17691453623978*rdxCp2Sq[0]*phiUx[1]+31.17691453623978*rdxCp2Sq[0]*phiLx[1]+(27.0*phiUx[0]+27.0*phiLx[0])*rdxCp2Sq[0])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = ((5.196152422706631*rdxCp2[1]*rho[3]+(30.0*rdxCp2[1]+60.0*rdxCp2[0])*rho[1])*volFac-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUx[3]+(69.28203230275508*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*phiLy[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLx[3]+90.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-90.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+60.0*phiLy[1]*rdxCp2Sq[1]+((-300.0*rdxCp2[0]*phiUx[1])+240.0*rdxCp2[0]*phiLy[1]-300.0*rdxCp2[0]*phiLx[1]+(259.8076211353315*phiUx[0]-259.8076211353315*phiLx[0])*rdxCp2[0])*rdxCp2[1]-600.0*rdxCp2Sq[0]*phiUx[1]-600.0*rdxCp2Sq[0]*phiLx[1]+(519.6152422706631*phiUx[0]-519.6152422706631*phiLx[0])*rdxCp2Sq[0])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 
  phiC[2] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[2]+17.32050807568877*rho[0]*rdxCp2[1])*volFac+((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]-60.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(25.98076211353316*phiUx[0]-51.96152422706631*phiLy[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((3.0*rdxCp2[1]+30.0*rdxCp2[0])*rho[3]+8.660254037844386*rdxCp2[1]*rho[1])*volFac+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+((-86.60254037844386*rdxCp2[0]*phiUx[1])-346.4101615137754*rdxCp2[0]*phiLy[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiUx[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(727.4613391789284*rdxCp2Sq[1]+4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+(4382.08854314926*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-43647.6803507357*rdxCp2R3[1])-269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(2182.384017536785*rdxCp2[0]*rdxCp2Sq[1]+16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(201600.0*rdxCp2R3[1]+1259760.0*rdxCp2[0]*rdxCp2Sq[1]+252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+37800.0*phiUy[0]*rdxCp2R3[1]+(16211.99555884469*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiUx[1]+(233370.0*phiUy[0]+31500.0*phiUx[0]+252000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2182.384017536785*rdxCp2Sq[0]*phiUy[1]-269472.4646415659*rdxCp2Sq[0]*phiUx[1]+(31500.0*phiUy[0]+233370.0*phiUx[0]+1259760.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-43647.6803507357*rdxCp2R3[0]*phiUx[1]+(37800.0*phiUx[0]+201600.0*bcVals[0])*rdxCp2R3[0])/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = (((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]+7274.613391789284*rdxCp2Sq[0]*rho[0])*volFac+((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-93600.0*rdxCp2[0]*rdxCp2Sq[1])-12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(12600.0*rdxCp2[0]*rdxCp2Sq[1]+8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(403221.4280020346*rdxCp2[0]*rdxCp2Sq[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+75600.0*phiUy[1]*rdxCp2R3[1]+(82980.0*rdxCp2[0]*phiUy[1]-210000.0*rdxCp2[0]*phiUx[1]+(81059.97779422344*phiUy[0]+181865.3347947321*phiUx[0]-1454922.678357857*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(11340.0*rdxCp2Sq[0]*phiUy[1]-341400.0*rdxCp2Sq[0]*phiUx[1]+(10911.92008768392*phiUy[0]+295661.0728520073*phiUx[0]-1313587.332460236*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-50400.0*rdxCp2R3[0]*phiUx[1]+(43647.6803507357*phiUx[0]-174590.7214029428*bcVals[0])*rdxCp2R3[0])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = (((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+7274.613391789284*rho[0]*rdxCp2Sq[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiUx[3]+((-50400.0*rdxCp2R3[1])-341400.0*rdxCp2[0]*rdxCp2Sq[1]-210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(11340.0*rdxCp2[0]*rdxCp2Sq[1]+82980.0*rdxCp2Sq[0]*rdxCp2[1]+75600.0*rdxCp2R3[0])*phiUx[2]+((-174590.7214029428*rdxCp2R3[1])-1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]-1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+43647.6803507357*phiUy[0]*rdxCp2R3[1]+(8640.0*rdxCp2[0]*phiUy[1]-12600.0*rdxCp2[0]*phiUx[1]+(295661.0728520073*phiUy[0]+10911.92008768392*phiUx[0]+87295.36070147139*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(12600.0*rdxCp2Sq[0]*phiUy[1]-93600.0*rdxCp2Sq[0]*phiUx[1]+(181865.3347947321*phiUy[0]+81059.97779422344*phiUx[0]+403221.4280020346*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(640.8587988004846*rdxCp2[0]*rdxCp2[1]+2424.871130596428*rdxCp2Sq[0])*rho[2]+(2424.871130596428*rdxCp2Sq[1]+640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]+5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiUx[3]+((-16627.68775266122*rdxCp2[0]*rdxCp2Sq[1])-24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(21823.84017536785*rdxCp2[0]*rdxCp2Sq[1]+128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]+29098.45356715714*rdxCp2R3[0])*phiUx[2]+(26400.0*rdxCp2[0]*rdxCp2Sq[1]-168000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+29098.45356715714*phiUy[1]*rdxCp2R3[1]+(128933.8621154272*rdxCp2[0]*phiUy[1]-24248.71130596428*rdxCp2[0]*phiUx[1]+(14400.0*phiUy[0]+21000.0*phiUx[0]-168000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]-16627.68775266122*rdxCp2Sq[0]*phiUx[1]+(21000.0*phiUy[0]+14400.0*phiUx[0]+26400.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(207.8460969082653*rdxCp2R3[1]+6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]+13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1143.153532995459*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-1091.192008768392*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*volFac+(2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(1385.640646055102*rdxCp2R4[1]+42816.29596310264*rdxCp2[0]*rdxCp2R3[1]+122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]+72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(623.5382907247956*rdxCp2[0]*rdxCp2R3[1]+22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]+48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(3200.0*rdxCp2R4[1]+96720.0*rdxCp2[0]*rdxCp2R3[1]+222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-1200.0*phiUy[0]*rdxCp2R4[1]+((-2078.460969082652*rdxCp2[0]*phiUy[1])+2078.460969082652*rdxCp2[0]*phiUx[1]+((-37080.0*phiUy[0])-1800.0*phiUx[0]-14400.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-6131.459858793824*rdxCp2Sq[0]*phiUy[1])+74720.67183852136*rdxCp2Sq[0]*phiUx[1]+((-106140.0*phiUy[0])-64710.0*phiUx[0]-359280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-4364.768035073569*rdxCp2R3[0]*phiUy[1])+188308.5637988883*rdxCp2R3[0]*phiUx[1]+((-63000.0*phiUy[0])-163080.0*phiUx[0]-879840.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+65471.52052610354*rdxCp2R4[0]*phiUx[1]+((-56700.0*phiUx[0])-302400.0*bcVals[0])*rdxCp2R4[0]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]-11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-10911.92008768392*rdxCp2R3[0]*rho[0])*volFac+(2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12000.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25200.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3600.0*rdxCp2[0]*rdxCp2R3[1]+25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(31869.73485926734*rdxCp2[0]*rdxCp2R3[1]+81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]+14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-2400.0*phiUy[1]*rdxCp2R4[1]+((-24960.0*rdxCp2[0]*phiUy[1])+12000.0*rdxCp2[0]*phiUx[1]+((-10392.30484541326*phiUy[0])-10392.30484541326*phiUx[0]+83138.4387633061*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-67140.0*rdxCp2Sq[0]*phiUy[1])+114600.0*rdxCp2Sq[0]*phiUx[1]+((-30657.29929396912*phiUy[0])-99246.51127369665*phiUx[0]+519615.2422706631*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-22680.0*rdxCp2R3[0]*phiUy[1])+237600.0*rdxCp2R3[0]*phiUx[1]+((-21823.84017536785*phiUy[0])-205767.6359391825*phiUx[0]+910365.9044582016*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+75600.0*rdxCp2R4[0]*phiUx[1]+(261886.0821044141*bcVals[0]-65471.52052610354*phiUx[0])*rdxCp2R4[0]))/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = (((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+(120.0*rdxCp2R3[1]+3870.0*rdxCp2[0]*rdxCp2Sq[1]+15150.0*rdxCp2Sq[0]*rdxCp2[1]+15750.0*rdxCp2R3[0])*rho[2]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiUx[3]+((-1800.0*rdxCp2[0]*rdxCp2R3[1])-50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-105000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(360.0*rdxCp2[0]*rdxCp2R3[1]+12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]+50760.0*rdxCp2R3[0]*rdxCp2[1]+56700.0*rdxCp2R4[0])*phiUx[2]+(1385.640646055102*rdxCp2R4[1]+43647.6803507357*rdxCp2[0]*rdxCp2R3[1]+145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]+121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-600.0*rdxCp2[0]*phiUy[1])+600.0*rdxCp2[0]*phiUx[1]+(1558.845726811989*phiUy[0]-519.6152422706631*phiUx[0]-4156.921938165305*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(21600.0*rdxCp2Sq[0]*phiUx[1]+(43647.6803507357*phiUy[0]-18706.14872174387*phiUx[0]-99766.1265159673*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(6300.0*rdxCp2R3[0]*phiUy[1]+46800.0*rdxCp2R3[0]*phiUx[1]+(90932.66739736605*phiUy[0]-40529.98889711172*phiUx[0]-201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiUx[3]+(3464.101615137754*rdxCp2[0]*rdxCp2R3[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(2078.460969082652*rdxCp2[0]*rdxCp2R3[1]+39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]+145907.9600296021*rdxCp2R3[0]*rdxCp2[1]+65471.52052610354*rdxCp2R4[0])*phiUx[2]+(10400.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(17320.50807568877*rdxCp2[0]*phiUy[1]+3464.101615137754*rdxCp2[0]*phiUx[1]+((-3000.0*phiUy[0])-3000.0*phiUx[0]+24000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(87295.36070147139*rdxCp2Sq[0]*phiUy[1]+24941.53162899183*rdxCp2Sq[0]*phiUx[1]+(86400.0*bcVals[0]-21600.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(32735.76026305177*rdxCp2R3[0]*phiUy[1]+24941.53162899183*rdxCp2R3[0]*phiUx[1]+(31500.0*phiUy[0]-21600.0*phiUx[0]-39600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1091.192008768392*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(13146.26562944778*rdxCp2[0]*rdxCp2Sq[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*volFac+((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(65471.52052610354*rdxCp2R4[1]+188308.5637988883*rdxCp2[0]*rdxCp2R3[1]+74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-4364.768035073569*rdxCp2[0]*rdxCp2R3[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-302400.0*rdxCp2R4[1])-879840.0*rdxCp2[0]*rdxCp2R3[1]-359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-56700.0*phiUy[0]*rdxCp2R4[1]+(48635.98667653406*rdxCp2[0]*phiUy[1]+72746.13391789283*rdxCp2[0]*phiUx[1]+((-163080.0*phiUy[0])-63000.0*phiUx[0]+42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(22447.37846609264*rdxCp2Sq[0]*phiUy[1]+122559.9151435737*rdxCp2Sq[0]*phiUx[1]+((-64710.0*phiUy[0])-106140.0*phiUx[0]+222560.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(623.5382907247956*rdxCp2R3[0]*phiUy[1]+42816.29596310264*rdxCp2R3[0]*phiUx[1]+((-1800.0*phiUy[0])-37080.0*phiUx[0]+96720.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+1385.640646055102*rdxCp2R4[0]*phiUx[1]+(3200.0*bcVals[0]-1200.0*phiUx[0])*rdxCp2R4[0]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = (((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-885.0*rdxCp2[0]*rdxCp2Sq[1])-360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(15750.0*rdxCp2R3[1]+15150.0*rdxCp2[0]*rdxCp2Sq[1]+3870.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(46800.0*rdxCp2[0]*rdxCp2R3[1]+21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(6300.0*rdxCp2[0]*rdxCp2R3[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+56700.0*phiUy[1]*rdxCp2R4[1]+(50760.0*rdxCp2[0]*phiUy[1]-105000.0*rdxCp2[0]*phiUx[1]+((-40529.98889711172*phiUy[0])+90932.66739736605*phiUx[0]+121243.5565298214*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(12870.0*rdxCp2Sq[0]*phiUy[1]-50400.0*rdxCp2Sq[0]*phiUx[1]+((-18706.14872174387*phiUy[0])+43647.6803507357*phiUx[0]+145838.6779972995*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(360.0*rdxCp2R3[0]*phiUy[1]-1800.0*rdxCp2R3[0]*phiUx[1]+((-519.6152422706631*phiUy[0])+1558.845726811989*phiUx[0]+43647.6803507357*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+1385.640646055102*bcVals[0]*rdxCp2R4[0])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-10911.92008768392*rho[0]*rdxCp2R3[1]-30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiUx[3]+(75600.0*rdxCp2R4[1]+237600.0*rdxCp2[0]*rdxCp2R3[1]+114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+12000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-22680.0*rdxCp2[0]*rdxCp2R3[1])-67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiUx[2]+(261886.0821044141*rdxCp2R4[1]+910365.9044582016*rdxCp2[0]*rdxCp2R3[1]+519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]+83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-65471.52052610354*phiUy[0]*rdxCp2R4[1]+(25920.0*rdxCp2[0]*phiUy[1]+25200.0*rdxCp2[0]*phiUx[1]+((-205767.6359391825*phiUy[0])-21823.84017536785*phiUx[0]+14549.22678357857*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(25920.0*rdxCp2Sq[0]*phiUy[1]+35400.0*rdxCp2Sq[0]*phiUx[1]+((-99246.51127369665*phiUy[0])-30657.29929396912*phiUx[0]+81752.798117251*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3600.0*rdxCp2R3[0]*phiUy[1]+12000.0*rdxCp2R3[0]*phiUx[1]+((-10392.30484541326*phiUy[0])-10392.30484541326*phiUx[0]+31869.73485926734*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(24941.53162899183*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(32735.76026305177*rdxCp2[0]*rdxCp2R3[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]+17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-39600.0*rdxCp2[0]*rdxCp2R3[1])+86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+65471.52052610354*phiUy[1]*rdxCp2R4[1]+(145907.9600296021*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiUx[1]+((-21600.0*phiUy[0])+31500.0*phiUx[0]+42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(39386.83536411626*rdxCp2Sq[0]*phiUy[1]+(35400.0*bcVals[0]-21600.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(2078.460969082652*rdxCp2R3[0]*phiUy[1]+3464.101615137754*rdxCp2R3[0]*phiUx[1]+((-3000.0*phiUy[0])-3000.0*phiUx[0]+10400.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-25.98076211353316*rdxCp2Sq[1])-285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+((-285.7883832488647*rdxCp2[0]*rdxCp2[1])-25.98076211353316*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*volFac+(600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-173.2050807568877*rdxCp2R3[1])-1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-103.9230484541326*rdxCp2[0]*rdxCp2Sq[1])-519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-400.0*rdxCp2R3[1])-4440.0*rdxCp2[0]*rdxCp2Sq[1]-200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+150.0*phiUy[0]*rdxCp2R3[1]+((-519.6152422706631*rdxCp2[0]*phiUy[1])-346.4101615137754*rdxCp2[0]*phiUx[1]+(1710.0*phiUy[0]+300.0*phiUx[0]-200.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2Sq[0]*phiUy[1])-1974.53792062852*rdxCp2Sq[0]*phiUx[1]+(300.0*phiUy[0]+1710.0*phiUx[0]-4440.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-173.2050807568877*rdxCp2R3[0]*phiUx[1]+(150.0*phiUx[0]-400.0*bcVals[0])*rdxCp2R3[0])/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0]*rho[0])*volFac+(1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(900.0*rdxCp2[0]*rdxCp2Sq[1]-900.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-7967.433714816835*rdxCp2[0]*rdxCp2Sq[1])-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-900.0*phiUy[1]*rdxCp2R3[1]+((-3060.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiUx[1]+(2598.076211353316*phiUy[0]-2598.076211353316*phiUx[0]-3464.101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-360.0*rdxCp2Sq[0]*phiUy[1])+600.0*rdxCp2Sq[0]*phiUx[1]+(519.6152422706631*phiUy[0]-519.6152422706631*phiUx[0]-12124.35565298214*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-1039.230484541326*bcVals[0]*rdxCp2R3[0]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+259.8076211353315*rho[0]*rdxCp2Sq[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiUx[3]+(600.0*rdxCp2[0]*rdxCp2Sq[1]+3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-3060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiUx[2]+((-1039.230484541326*rdxCp2R3[1])-12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-900.0*rdxCp2[0]*phiUy[1])-600.0*rdxCp2[0]*phiUx[1]+((-519.6152422706631*phiUy[0])+519.6152422706631*phiUx[0]-346.4101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(900.0*rdxCp2Sq[0]*phiUy[1]-3000.0*rdxCp2Sq[0]*phiUx[1]+((-2598.076211353316*phiUy[0])+2598.076211353316*phiUx[0]-7967.433714816835*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-129.9038105676658*rdxCp2Sq[0])*rho[2]+((-129.9038105676658*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(866.0254037844386*rdxCp2Sq[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(519.6152422706631*rdxCp2[0]*rdxCp2Sq[1]+2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-2600.0*rdxCp2[0]*rdxCp2Sq[1])-1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(2598.076211353316*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiUx[1]+(750.0*phiUy[0]-750.0*phiUx[0]-1000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2Sq[0]*phiUy[1]-866.0254037844386*rdxCp2Sq[0]*phiUx[1]+((-750.0*phiUy[0])+750.0*phiUx[0]-2600.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(727.4613391789284*rdxCp2Sq[1]+4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4382.08854314926*rdxCp2[0]*rdxCp2[1])-727.4613391789284*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-201600.0*rdxCp2R3[1])-1259760.0*rdxCp2[0]*rdxCp2Sq[1]-252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2182.384017536785*rdxCp2[0]*rdxCp2Sq[1]+16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-43647.6803507357*rdxCp2R3[1])-269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-37800.0*phiLy[0]*rdxCp2R3[1]+(36373.06695894642*rdxCp2[0]*phiUx[1]-16211.99555884469*rdxCp2[0]*phiLy[1]+((-31500.0*phiUx[0])-233370.0*phiLy[0]-252000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(269472.4646415659*rdxCp2Sq[0]*phiUx[1]-2182.384017536785*rdxCp2Sq[0]*phiLy[1]+((-233370.0*phiUx[0])-31500.0*phiLy[0]-1259760.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+43647.6803507357*rdxCp2R3[0]*phiUx[1]+((-37800.0*phiUx[0])-201600.0*bcVals[0])*rdxCp2R3[0]))/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]-7274.613391789284*rdxCp2Sq[0]*rho[0])*volFac+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-403221.4280020346*rdxCp2[0]*rdxCp2Sq[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(12600.0*rdxCp2[0]*rdxCp2Sq[1]+8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-93600.0*rdxCp2[0]*rdxCp2Sq[1])-12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-75600.0*phiLy[1]*rdxCp2R3[1]+(210000.0*rdxCp2[0]*phiUx[1]-82980.0*rdxCp2[0]*phiLy[1]+((-181865.3347947321*phiUx[0])-81059.97779422344*phiLy[0]+1454922.678357857*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(341400.0*rdxCp2Sq[0]*phiUx[1]-11340.0*rdxCp2Sq[0]*phiLy[1]+((-295661.0728520073*phiUx[0])-10911.92008768392*phiLy[0]+1313587.332460236*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+50400.0*rdxCp2R3[0]*phiUx[1]+(174590.7214029428*bcVals[0]-43647.6803507357*phiUx[0])*rdxCp2R3[0]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = (((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-7274.613391789284*rho[0]*rdxCp2Sq[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiUx[3]+((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(174590.7214029428*rdxCp2R3[1]+1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]+1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(11340.0*rdxCp2[0]*rdxCp2Sq[1]+82980.0*rdxCp2Sq[0]*rdxCp2[1]+75600.0*rdxCp2R3[0])*phiUx[2]+((-50400.0*rdxCp2R3[1])-341400.0*rdxCp2[0]*rdxCp2Sq[1]-210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-43647.6803507357*phiLy[0]*rdxCp2R3[1]+(12600.0*rdxCp2[0]*phiUx[1]-8640.0*rdxCp2[0]*phiLy[1]+((-10911.92008768392*phiUx[0])-295661.0728520073*phiLy[0]-87295.36070147139*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(93600.0*rdxCp2Sq[0]*phiUx[1]-12600.0*rdxCp2Sq[0]*phiLy[1]+((-81059.97779422344*phiUx[0])-181865.3347947321*phiLy[0]-403221.4280020346*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(640.8587988004846*rdxCp2[0]*rdxCp2[1]+2424.871130596428*rdxCp2Sq[0])*rho[2]+((-2424.871130596428*rdxCp2Sq[1])-640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]-5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiUx[3]+((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(168000.0*rdxCp2Sq[0]*rdxCp2[1]-26400.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(21823.84017536785*rdxCp2[0]*rdxCp2Sq[1]+128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]+29098.45356715714*rdxCp2R3[0])*phiUx[2]+((-16627.68775266122*rdxCp2[0]*rdxCp2Sq[1])-24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-29098.45356715714*phiLy[1]*rdxCp2R3[1]+(24248.71130596428*rdxCp2[0]*phiUx[1]-128933.8621154272*rdxCp2[0]*phiLy[1]+((-21000.0*phiUx[0])-14400.0*phiLy[0]+168000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(16627.68775266122*rdxCp2Sq[0]*phiUx[1]-21823.84017536785*rdxCp2Sq[0]*phiLy[1]+((-14400.0*phiUx[0])-21000.0*phiLy[0]-26400.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(207.8460969082653*rdxCp2R3[1]+6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]+13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1143.153532995459*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+1091.192008768392*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*volFac+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3200.0*rdxCp2R4[1]+96720.0*rdxCp2[0]*rdxCp2R3[1]+222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(623.5382907247956*rdxCp2[0]*rdxCp2R3[1]+22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]+48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(1385.640646055102*rdxCp2R4[1]+42816.29596310264*rdxCp2[0]*rdxCp2R3[1]+122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]+72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+1200.0*phiLy[0]*rdxCp2R4[1]+((-2078.460969082652*rdxCp2[0]*phiUx[1])+2078.460969082652*rdxCp2[0]*phiLy[1]+(1800.0*phiUx[0]+37080.0*phiLy[0]+14400.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-74720.67183852136*rdxCp2Sq[0]*phiUx[1])+6131.459858793824*rdxCp2Sq[0]*phiLy[1]+(64710.0*phiUx[0]+106140.0*phiLy[0]+359280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-188308.5637988883*rdxCp2R3[0]*phiUx[1])+4364.768035073569*rdxCp2R3[0]*phiLy[1]+(163080.0*phiUx[0]+63000.0*phiLy[0]+879840.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-65471.52052610354*rdxCp2R4[0]*phiUx[1]+(56700.0*phiUx[0]+302400.0*bcVals[0])*rdxCp2R4[0])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = (((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]+11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+10911.92008768392*rdxCp2R3[0]*rho[0])*volFac+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(31869.73485926734*rdxCp2[0]*rdxCp2R3[1]+81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]+14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3600.0*rdxCp2[0]*rdxCp2R3[1]+25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(12000.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25200.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+2400.0*phiLy[1]*rdxCp2R4[1]+((-12000.0*rdxCp2[0]*phiUx[1])+24960.0*rdxCp2[0]*phiLy[1]+(10392.30484541326*phiUx[0]+10392.30484541326*phiLy[0]-83138.4387633061*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-114600.0*rdxCp2Sq[0]*phiUx[1])+67140.0*rdxCp2Sq[0]*phiLy[1]+(99246.51127369665*phiUx[0]+30657.29929396912*phiLy[0]-519615.2422706631*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-237600.0*rdxCp2R3[0]*phiUx[1])+22680.0*rdxCp2R3[0]*phiLy[1]+(205767.6359391825*phiUx[0]+21823.84017536785*phiLy[0]-910365.9044582016*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-75600.0*rdxCp2R4[0]*phiUx[1]+(65471.52052610354*phiUx[0]-261886.0821044141*bcVals[0])*rdxCp2R4[0])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = (((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+(120.0*rdxCp2R3[1]+3870.0*rdxCp2[0]*rdxCp2Sq[1]+15150.0*rdxCp2Sq[0]*rdxCp2[1]+15750.0*rdxCp2R3[0])*rho[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiUx[3]+(692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(1385.640646055102*rdxCp2R4[1]+43647.6803507357*rdxCp2[0]*rdxCp2R3[1]+145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]+121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(360.0*rdxCp2[0]*rdxCp2R3[1]+12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]+50760.0*rdxCp2R3[0]*rdxCp2[1]+56700.0*rdxCp2R4[0])*phiUx[2]+((-1800.0*rdxCp2[0]*rdxCp2R3[1])-50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-105000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-600.0*rdxCp2[0]*phiUx[1])+600.0*rdxCp2[0]*phiLy[1]+(519.6152422706631*phiUx[0]-1558.845726811989*phiLy[0]+4156.921938165305*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((18706.14872174387*phiUx[0]-43647.6803507357*phiLy[0]+99766.1265159673*bcVals[0])*rdxCp2Sq[0]-21600.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-46800.0*rdxCp2R3[0]*phiUx[1])-6300.0*rdxCp2R3[0]*phiLy[1]+(40529.98889711172*phiUx[0]-90932.66739736605*phiLy[0]+201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiUx[3]+((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(10400.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(2078.460969082652*rdxCp2[0]*rdxCp2R3[1]+39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]+145907.9600296021*rdxCp2R3[0]*rdxCp2[1]+65471.52052610354*rdxCp2R4[0])*phiUx[2]+(3464.101615137754*rdxCp2[0]*rdxCp2R3[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-3464.101615137754*rdxCp2[0]*phiUx[1])-17320.50807568877*rdxCp2[0]*phiLy[1]+(3000.0*phiUx[0]+3000.0*phiLy[0]-24000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-24941.53162899183*rdxCp2Sq[0]*phiUx[1])-87295.36070147139*rdxCp2Sq[0]*phiLy[1]+(21600.0*phiUx[0]-86400.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-24941.53162899183*rdxCp2R3[0]*phiUx[1])-32735.76026305177*rdxCp2R3[0]*phiLy[1]+(21600.0*phiUx[0]-31500.0*phiLy[0]+39600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1091.192008768392*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-13146.26562944778*rdxCp2[0]*rdxCp2Sq[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]-207.8460969082653*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*volFac+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(302400.0*rdxCp2R4[1]+879840.0*rdxCp2[0]*rdxCp2R3[1]+359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-4364.768035073569*rdxCp2[0]*rdxCp2R3[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(65471.52052610354*rdxCp2R4[1]+188308.5637988883*rdxCp2[0]*rdxCp2R3[1]+74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+56700.0*phiLy[0]*rdxCp2R4[1]+((-72746.13391789283*rdxCp2[0]*phiUx[1])-48635.98667653406*rdxCp2[0]*phiLy[1]+(63000.0*phiUx[0]+163080.0*phiLy[0]-42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-122559.9151435737*rdxCp2Sq[0]*phiUx[1])-22447.37846609264*rdxCp2Sq[0]*phiLy[1]+(106140.0*phiUx[0]+64710.0*phiLy[0]-222560.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-42816.29596310264*rdxCp2R3[0]*phiUx[1])-623.5382907247956*rdxCp2R3[0]*phiLy[1]+(37080.0*phiUx[0]+1800.0*phiLy[0]-96720.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-1385.640646055102*rdxCp2R4[0]*phiUx[1]+(1200.0*phiUx[0]-3200.0*bcVals[0])*rdxCp2R4[0])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-885.0*rdxCp2[0]*rdxCp2Sq[1])-360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-15750.0*rdxCp2R3[1])-15150.0*rdxCp2[0]*rdxCp2Sq[1]-3870.0*rdxCp2Sq[0]*rdxCp2[1]-120.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]+4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(6300.0*rdxCp2[0]*rdxCp2R3[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(46800.0*rdxCp2[0]*rdxCp2R3[1]+21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-56700.0*phiLy[1]*rdxCp2R4[1]+(105000.0*rdxCp2[0]*phiUx[1]-50760.0*rdxCp2[0]*phiLy[1]+((-90932.66739736605*phiUx[0])+40529.98889711172*phiLy[0]-121243.5565298214*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(50400.0*rdxCp2Sq[0]*phiUx[1]-12870.0*rdxCp2Sq[0]*phiLy[1]+((-43647.6803507357*phiUx[0])+18706.14872174387*phiLy[0]-145838.6779972995*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(1800.0*rdxCp2R3[0]*phiUx[1]-360.0*rdxCp2R3[0]*phiLy[1]+((-1558.845726811989*phiUx[0])+519.6152422706631*phiLy[0]-43647.6803507357*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-1385.640646055102*bcVals[0]*rdxCp2R4[0]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+10911.92008768392*rho[0]*rdxCp2R3[1]+30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiUx[3]+((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-261886.0821044141*rdxCp2R4[1])-910365.9044582016*rdxCp2[0]*rdxCp2R3[1]-519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]-83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-22680.0*rdxCp2[0]*rdxCp2R3[1])-67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiUx[2]+(75600.0*rdxCp2R4[1]+237600.0*rdxCp2[0]*rdxCp2R3[1]+114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+12000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+65471.52052610354*phiLy[0]*rdxCp2R4[1]+((-25200.0*rdxCp2[0]*phiUx[1])-25920.0*rdxCp2[0]*phiLy[1]+(21823.84017536785*phiUx[0]+205767.6359391825*phiLy[0]-14549.22678357857*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-35400.0*rdxCp2Sq[0]*phiUx[1])-25920.0*rdxCp2Sq[0]*phiLy[1]+(30657.29929396912*phiUx[0]+99246.51127369665*phiLy[0]-81752.798117251*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-12000.0*rdxCp2R3[0]*phiUx[1])-3600.0*rdxCp2R3[0]*phiLy[1]+(10392.30484541326*phiUx[0]+10392.30484541326*phiLy[0]-31869.73485926734*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]))/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(39600.0*rdxCp2[0]*rdxCp2R3[1]-86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(32735.76026305177*rdxCp2[0]*rdxCp2R3[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]+17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(24941.53162899183*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-65471.52052610354*phiLy[1]*rdxCp2R4[1]+(36373.06695894642*rdxCp2[0]*phiUx[1]-145907.9600296021*rdxCp2[0]*phiLy[1]+((-31500.0*phiUx[0])+21600.0*phiLy[0]-42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((21600.0*phiLy[0]-35400.0*bcVals[0])*rdxCp2Sq[0]-39386.83536411626*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-3464.101615137754*rdxCp2R3[0]*phiUx[1])-2078.460969082652*rdxCp2R3[0]*phiLy[1]+(3000.0*phiUx[0]+3000.0*phiLy[0]-10400.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-25.98076211353316*rdxCp2Sq[1])-285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+(285.7883832488647*rdxCp2[0]*rdxCp2[1]+25.98076211353316*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*volFac+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-400.0*rdxCp2R3[1])-4440.0*rdxCp2[0]*rdxCp2Sq[1]-200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-103.9230484541326*rdxCp2[0]*rdxCp2Sq[1])-519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-173.2050807568877*rdxCp2R3[1])-1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-150.0*phiLy[0]*rdxCp2R3[1]+(346.4101615137754*rdxCp2[0]*phiUx[1]+519.6152422706631*rdxCp2[0]*phiLy[1]+((-300.0*phiUx[0])-1710.0*phiLy[0]+200.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1974.53792062852*rdxCp2Sq[0]*phiUx[1]+103.9230484541326*rdxCp2Sq[0]*phiLy[1]+((-1710.0*phiUx[0])-300.0*phiLy[0]+4440.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+173.2050807568877*rdxCp2R3[0]*phiUx[1]+(400.0*bcVals[0]-150.0*phiUx[0])*rdxCp2R3[0]))/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = (((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]-259.8076211353315*rdxCp2Sq[0]*rho[0])*volFac+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-7967.433714816835*rdxCp2[0]*rdxCp2Sq[1])-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(900.0*rdxCp2[0]*rdxCp2Sq[1]-900.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+900.0*phiLy[1]*rdxCp2R3[1]+((-3000.0*rdxCp2[0]*phiUx[1])+3060.0*rdxCp2[0]*phiLy[1]+(2598.076211353316*phiUx[0]-2598.076211353316*phiLy[0]+3464.101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-600.0*rdxCp2Sq[0]*phiUx[1])+360.0*rdxCp2Sq[0]*phiLy[1]+(519.6152422706631*phiUx[0]-519.6152422706631*phiLy[0]+12124.35565298214*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+1039.230484541326*bcVals[0]*rdxCp2R3[0])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-259.8076211353315*rho[0]*rdxCp2Sq[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiUx[3]+(1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1039.230484541326*rdxCp2R3[1])-12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-3060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiUx[2]+(600.0*rdxCp2[0]*rdxCp2Sq[1]+3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(600.0*rdxCp2[0]*phiUx[1]+900.0*rdxCp2[0]*phiLy[1]+((-519.6152422706631*phiUx[0])+519.6152422706631*phiLy[0]+346.4101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(3000.0*rdxCp2Sq[0]*phiUx[1]-900.0*rdxCp2Sq[0]*phiLy[1]+((-2598.076211353316*phiUx[0])+2598.076211353316*phiLy[0]+7967.433714816835*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-129.9038105676658*rdxCp2Sq[0])*rho[2]+(129.9038105676658*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-2600.0*rdxCp2[0]*rdxCp2Sq[1])-1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(519.6152422706631*rdxCp2[0]*rdxCp2Sq[1]+2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(866.0254037844386*rdxCp2Sq[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-866.0254037844386*rdxCp2[0]*phiUx[1])-2598.076211353316*rdxCp2[0]*phiLy[1]+(750.0*phiUx[0]-750.0*phiLy[0]+1000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(866.0254037844386*rdxCp2Sq[0]*phiUx[1]-519.6152422706631*rdxCp2Sq[0]*phiLy[1]+((-750.0*phiUx[0])+750.0*phiLy[0]+2600.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-727.4613391789284*rdxCp2Sq[1])-4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+(4382.08854314926*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(43647.6803507357*rdxCp2R3[1]+269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-2182.384017536785*rdxCp2[0]*rdxCp2Sq[1])-16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-201600.0*rdxCp2R3[1])-1259760.0*rdxCp2[0]*rdxCp2Sq[1]-252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-37800.0*phiUy[0]*rdxCp2R3[1]+(16211.99555884469*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiLx[1]-252000.0*rdxCp2[0]*bcVals[1]+((-233370.0*phiUy[0])-31500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(2182.384017536785*rdxCp2Sq[0]*phiUy[1]-269472.4646415659*rdxCp2Sq[0]*phiLx[1]-1259760.0*rdxCp2Sq[0]*bcVals[1]+((-31500.0*phiUy[0])-233370.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-43647.6803507357*rdxCp2R3[0]*phiLx[1]-201600.0*rdxCp2R3[0]*bcVals[1]-37800.0*phiLx[0]*rdxCp2R3[0]))/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = (((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]-7274.613391789284*rdxCp2Sq[0]*rho[0])*volFac+((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(93600.0*rdxCp2[0]*rdxCp2Sq[1]+12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-12600.0*rdxCp2[0]*rdxCp2Sq[1])-8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-403221.4280020346*rdxCp2[0]*rdxCp2Sq[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+75600.0*phiUy[1]*rdxCp2R3[1]+(82980.0*rdxCp2[0]*phiUy[1]-210000.0*rdxCp2[0]*phiLx[1]+1454922.678357857*rdxCp2[0]*bcVals[1]+((-81059.97779422344*phiUy[0])-181865.3347947321*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(11340.0*rdxCp2Sq[0]*phiUy[1]-341400.0*rdxCp2Sq[0]*phiLx[1]+1313587.332460236*rdxCp2Sq[0]*bcVals[1]+((-10911.92008768392*phiUy[0])-295661.0728520073*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-50400.0*rdxCp2R3[0]*phiLx[1]+174590.7214029428*rdxCp2R3[0]*bcVals[1]-43647.6803507357*phiLx[0]*rdxCp2R3[0])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-7274.613391789284*rho[0]*rdxCp2Sq[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiLx[3]+(50400.0*rdxCp2R3[1]+341400.0*rdxCp2[0]*rdxCp2Sq[1]+210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-11340.0*rdxCp2[0]*rdxCp2Sq[1])-82980.0*rdxCp2Sq[0]*rdxCp2[1]-75600.0*rdxCp2R3[0])*phiLx[2]+(174590.7214029428*rdxCp2R3[1]+1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]+1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-43647.6803507357*phiUy[0]*rdxCp2R3[1]+(8640.0*rdxCp2[0]*phiUy[1]-12600.0*rdxCp2[0]*phiLx[1]-87295.36070147139*rdxCp2[0]*bcVals[1]+((-295661.0728520073*phiUy[0])-10911.92008768392*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(12600.0*rdxCp2Sq[0]*phiUy[1]-93600.0*rdxCp2Sq[0]*phiLx[1]-403221.4280020346*rdxCp2Sq[0]*bcVals[1]+((-181865.3347947321*phiUy[0])-81059.97779422344*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-640.8587988004846*rdxCp2[0]*rdxCp2[1])-2424.871130596428*rdxCp2Sq[0])*rho[2]+(2424.871130596428*rdxCp2Sq[1]+640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]-5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiLx[3]+(16627.68775266122*rdxCp2[0]*rdxCp2Sq[1]+24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-21823.84017536785*rdxCp2[0]*rdxCp2Sq[1])-128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]-29098.45356715714*rdxCp2R3[0])*phiLx[2]+(168000.0*rdxCp2Sq[0]*rdxCp2[1]-26400.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+29098.45356715714*phiUy[1]*rdxCp2R3[1]+(128933.8621154272*rdxCp2[0]*phiUy[1]-24248.71130596428*rdxCp2[0]*phiLx[1]+168000.0*rdxCp2[0]*bcVals[1]+((-14400.0*phiUy[0])-21000.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]-16627.68775266122*rdxCp2Sq[0]*phiLx[1]-26400.0*rdxCp2Sq[0]*bcVals[1]+((-21000.0*phiUy[0])-14400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-207.8460969082653*rdxCp2R3[1])-6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]-13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1143.153532995459*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-1091.192008768392*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*volFac+(2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-1385.640646055102*rdxCp2R4[1])-42816.29596310264*rdxCp2[0]*rdxCp2R3[1]-122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]-72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-623.5382907247956*rdxCp2[0]*rdxCp2R3[1])-22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]-48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-3200.0*rdxCp2R4[1])-96720.0*rdxCp2[0]*rdxCp2R3[1]-222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+1200.0*phiUy[0]*rdxCp2R4[1]+((-2078.460969082652*rdxCp2[0]*phiUy[1])+2078.460969082652*rdxCp2[0]*phiLx[1]+14400.0*rdxCp2[0]*bcVals[1]+(37080.0*phiUy[0]+1800.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-6131.459858793824*rdxCp2Sq[0]*phiUy[1])+74720.67183852136*rdxCp2Sq[0]*phiLx[1]+359280.0*rdxCp2Sq[0]*bcVals[1]+(106140.0*phiUy[0]+64710.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-4364.768035073569*rdxCp2R3[0]*phiUy[1])+188308.5637988883*rdxCp2R3[0]*phiLx[1]+879840.0*rdxCp2R3[0]*bcVals[1]+(63000.0*phiUy[0]+163080.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+65471.52052610354*rdxCp2R4[0]*phiLx[1]+302400.0*rdxCp2R4[0]*bcVals[1]+56700.0*phiLx[0]*rdxCp2R4[0])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]+11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+10911.92008768392*rdxCp2R3[0]*rho[0])*volFac+(2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-12000.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25200.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-3600.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-31869.73485926734*rdxCp2[0]*rdxCp2R3[1])-81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]-14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-2400.0*phiUy[1]*rdxCp2R4[1]+((-24960.0*rdxCp2[0]*phiUy[1])+12000.0*rdxCp2[0]*phiLx[1]-83138.4387633061*rdxCp2[0]*bcVals[1]+(10392.30484541326*phiUy[0]+10392.30484541326*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-67140.0*rdxCp2Sq[0]*phiUy[1])+114600.0*rdxCp2Sq[0]*phiLx[1]-519615.2422706631*rdxCp2Sq[0]*bcVals[1]+(30657.29929396912*phiUy[0]+99246.51127369665*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-22680.0*rdxCp2R3[0]*phiUy[1])+237600.0*rdxCp2R3[0]*phiLx[1]-910365.9044582016*rdxCp2R3[0]*bcVals[1]+(21823.84017536785*phiUy[0]+205767.6359391825*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+75600.0*rdxCp2R4[0]*phiLx[1]-261886.0821044141*rdxCp2R4[0]*bcVals[1]+65471.52052610354*phiLx[0]*rdxCp2R4[0]))/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+((-120.0*rdxCp2R3[1])-3870.0*rdxCp2[0]*rdxCp2Sq[1]-15150.0*rdxCp2Sq[0]*rdxCp2[1]-15750.0*rdxCp2R3[0])*rho[2]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[3]+(1800.0*rdxCp2[0]*rdxCp2R3[1]+50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+105000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]-50760.0*rdxCp2R3[0]*rdxCp2[1]-56700.0*rdxCp2R4[0])*phiLx[2]+((-1385.640646055102*rdxCp2R4[1])-43647.6803507357*rdxCp2[0]*rdxCp2R3[1]-145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]-121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-600.0*rdxCp2[0]*phiUy[1])+600.0*rdxCp2[0]*phiLx[1]+4156.921938165305*rdxCp2[0]*bcVals[1]+(519.6152422706631*phiLx[0]-1558.845726811989*phiUy[0])*rdxCp2[0])*rdxCp2R3[1]+(21600.0*rdxCp2Sq[0]*phiLx[1]+99766.1265159673*rdxCp2Sq[0]*bcVals[1]+(18706.14872174387*phiLx[0]-43647.6803507357*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(6300.0*rdxCp2R3[0]*phiUy[1]+46800.0*rdxCp2R3[0]*phiLx[1]+201610.7140010173*rdxCp2R3[0]*bcVals[1]+(40529.98889711172*phiLx[0]-90932.66739736605*phiUy[0])*rdxCp2R3[0])*rdxCp2[1]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiLx[3]+(36373.06695894642*rdxCp2R3[0]*rdxCp2[1]-3464.101615137754*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-2078.460969082652*rdxCp2[0]*rdxCp2R3[1])-39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]-145907.9600296021*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[2]+((-10400.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(17320.50807568877*rdxCp2[0]*phiUy[1]+3464.101615137754*rdxCp2[0]*phiLx[1]-24000.0*rdxCp2[0]*bcVals[1]+(3000.0*phiUy[0]+3000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(87295.36070147139*rdxCp2Sq[0]*phiUy[1]+24941.53162899183*rdxCp2Sq[0]*phiLx[1]-86400.0*rdxCp2Sq[0]*bcVals[1]+21600.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(32735.76026305177*rdxCp2R3[0]*phiUy[1]+24941.53162899183*rdxCp2R3[0]*phiLx[1]+39600.0*rdxCp2R3[0]*bcVals[1]+(21600.0*phiLx[0]-31500.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1091.192008768392*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(13146.26562944778*rdxCp2[0]*rdxCp2Sq[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*volFac+((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-65471.52052610354*rdxCp2R4[1])-188308.5637988883*rdxCp2[0]*rdxCp2R3[1]-74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(4364.768035073569*rdxCp2[0]*rdxCp2R3[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(302400.0*rdxCp2R4[1]+879840.0*rdxCp2[0]*rdxCp2R3[1]+359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+56700.0*phiUy[0]*rdxCp2R4[1]+(48635.98667653406*rdxCp2[0]*phiUy[1]+72746.13391789283*rdxCp2[0]*phiLx[1]+42000.0*rdxCp2[0]*bcVals[1]+(163080.0*phiUy[0]+63000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(22447.37846609264*rdxCp2Sq[0]*phiUy[1]+122559.9151435737*rdxCp2Sq[0]*phiLx[1]+222560.0*rdxCp2Sq[0]*bcVals[1]+(64710.0*phiUy[0]+106140.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(623.5382907247956*rdxCp2R3[0]*phiUy[1]+42816.29596310264*rdxCp2R3[0]*phiLx[1]+96720.0*rdxCp2R3[0]*bcVals[1]+(1800.0*phiUy[0]+37080.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+1385.640646055102*rdxCp2R4[0]*phiLx[1]+3200.0*rdxCp2R4[0]*bcVals[1]+1200.0*phiLx[0]*rdxCp2R4[0])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = (((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(885.0*rdxCp2[0]*rdxCp2Sq[1]+360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(15750.0*rdxCp2R3[1]+15150.0*rdxCp2[0]*rdxCp2Sq[1]+3870.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-46800.0*rdxCp2[0]*rdxCp2R3[1])-21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(600.0*rdxCp2R3[0]*rdxCp2[1]-6300.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]+4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+56700.0*phiUy[1]*rdxCp2R4[1]+(50760.0*rdxCp2[0]*phiUy[1]-105000.0*rdxCp2[0]*phiLx[1]+121243.5565298214*rdxCp2[0]*bcVals[1]+(40529.98889711172*phiUy[0]-90932.66739736605*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(12870.0*rdxCp2Sq[0]*phiUy[1]-50400.0*rdxCp2Sq[0]*phiLx[1]+145838.6779972995*rdxCp2Sq[0]*bcVals[1]+(18706.14872174387*phiUy[0]-43647.6803507357*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(360.0*rdxCp2R3[0]*phiUy[1]-1800.0*rdxCp2R3[0]*phiLx[1]+43647.6803507357*rdxCp2R3[0]*bcVals[1]+(519.6152422706631*phiUy[0]-1558.845726811989*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+1385.640646055102*rdxCp2R4[0]*bcVals[1])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = (((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+10911.92008768392*rho[0]*rdxCp2R3[1]+30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiLx[3]+((-75600.0*rdxCp2R4[1])-237600.0*rdxCp2[0]*rdxCp2R3[1]-114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-12000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(22680.0*rdxCp2[0]*rdxCp2R3[1]+67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiLx[2]+((-261886.0821044141*rdxCp2R4[1])-910365.9044582016*rdxCp2[0]*rdxCp2R3[1]-519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]-83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+65471.52052610354*phiUy[0]*rdxCp2R4[1]+(25920.0*rdxCp2[0]*phiUy[1]+25200.0*rdxCp2[0]*phiLx[1]+14549.22678357857*rdxCp2[0]*bcVals[1]+(205767.6359391825*phiUy[0]+21823.84017536785*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(25920.0*rdxCp2Sq[0]*phiUy[1]+35400.0*rdxCp2Sq[0]*phiLx[1]+81752.798117251*rdxCp2Sq[0]*bcVals[1]+(99246.51127369665*phiUy[0]+30657.29929396912*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3600.0*rdxCp2R3[0]*phiUy[1]+12000.0*rdxCp2R3[0]*phiLx[1]+31869.73485926734*rdxCp2R3[0]*bcVals[1]+(10392.30484541326*phiUy[0]+10392.30484541326*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-24941.53162899183*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-32735.76026305177*rdxCp2[0]*rdxCp2R3[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]-17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(39600.0*rdxCp2[0]*rdxCp2R3[1]-86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+65471.52052610354*phiUy[1]*rdxCp2R4[1]+(145907.9600296021*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiLx[1]+42000.0*rdxCp2[0]*bcVals[1]+(21600.0*phiUy[0]-31500.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(39386.83536411626*rdxCp2Sq[0]*phiUy[1]+35400.0*rdxCp2Sq[0]*bcVals[1]+21600.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(2078.460969082652*rdxCp2R3[0]*phiUy[1]+3464.101615137754*rdxCp2R3[0]*phiLx[1]+10400.0*rdxCp2R3[0]*bcVals[1]+(3000.0*phiUy[0]+3000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(25.98076211353316*rdxCp2Sq[1]+285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+((-285.7883832488647*rdxCp2[0]*rdxCp2[1])-25.98076211353316*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*volFac+(600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(173.2050807568877*rdxCp2R3[1]+1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(400.0*rdxCp2R3[1]+4440.0*rdxCp2[0]*rdxCp2Sq[1]+200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-150.0*phiUy[0]*rdxCp2R3[1]+((-519.6152422706631*rdxCp2[0]*phiUy[1])-346.4101615137754*rdxCp2[0]*phiLx[1]-200.0*rdxCp2[0]*bcVals[1]+((-1710.0*phiUy[0])-300.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2Sq[0]*phiUy[1])-1974.53792062852*rdxCp2Sq[0]*phiLx[1]-4440.0*rdxCp2Sq[0]*bcVals[1]+((-300.0*phiUy[0])-1710.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-173.2050807568877*rdxCp2R3[0]*phiLx[1]-400.0*rdxCp2R3[0]*bcVals[1]-150.0*phiLx[0]*rdxCp2R3[0]))/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]-259.8076211353315*rdxCp2Sq[0]*rho[0])*volFac+(1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(3000.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(900.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(7967.433714816835*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-900.0*phiUy[1]*rdxCp2R3[1]+((-3060.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiLx[1]-3464.101615137754*rdxCp2[0]*bcVals[1]+(2598.076211353316*phiLx[0]-2598.076211353316*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-360.0*rdxCp2Sq[0]*phiUy[1])+600.0*rdxCp2Sq[0]*phiLx[1]-12124.35565298214*rdxCp2Sq[0]*bcVals[1]+(519.6152422706631*phiLx[0]-519.6152422706631*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]-1039.230484541326*rdxCp2R3[0]*bcVals[1]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = (((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-259.8076211353315*rho[0]*rdxCp2Sq[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiLx[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+3060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiLx[2]+(1039.230484541326*rdxCp2R3[1]+12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-900.0*rdxCp2[0]*phiUy[1])-600.0*rdxCp2[0]*phiLx[1]-346.4101615137754*rdxCp2[0]*bcVals[1]+(519.6152422706631*phiUy[0]-519.6152422706631*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(900.0*rdxCp2Sq[0]*phiUy[1]-3000.0*rdxCp2Sq[0]*phiLx[1]-7967.433714816835*rdxCp2Sq[0]*bcVals[1]+(2598.076211353316*phiUy[0]-2598.076211353316*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0])*rho[2]+((-129.9038105676658*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(866.0254037844386*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2Sq[1])-2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(2600.0*rdxCp2[0]*rdxCp2Sq[1]+1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(2598.076211353316*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiLx[1]-1000.0*rdxCp2[0]*bcVals[1]+(750.0*phiLx[0]-750.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2Sq[0]*phiUy[1]-866.0254037844386*rdxCp2Sq[0]*phiLx[1]-2600.0*rdxCp2Sq[0]*bcVals[1]+(750.0*phiUy[0]-750.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-727.4613391789284*rdxCp2Sq[1])-4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4382.08854314926*rdxCp2[0]*rdxCp2[1])-727.4613391789284*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*volFac+((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(201600.0*rdxCp2R3[1]+1259760.0*rdxCp2[0]*rdxCp2Sq[1]+252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(43647.6803507357*rdxCp2R3[1]+269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-2182.384017536785*rdxCp2[0]*rdxCp2Sq[1])-16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+37800.0*phiLy[0]*rdxCp2R3[1]+((-16211.99555884469*rdxCp2[0]*phiLy[1])+36373.06695894642*rdxCp2[0]*phiLx[1]+252000.0*rdxCp2[0]*bcVals[1]+(233370.0*phiLy[0]+31500.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-2182.384017536785*rdxCp2Sq[0]*phiLy[1])+269472.4646415659*rdxCp2Sq[0]*phiLx[1]+1259760.0*rdxCp2Sq[0]*bcVals[1]+(31500.0*phiLy[0]+233370.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+43647.6803507357*rdxCp2R3[0]*phiLx[1]+201600.0*rdxCp2R3[0]*bcVals[1]+37800.0*phiLx[0]*rdxCp2R3[0])/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]+7274.613391789284*rdxCp2Sq[0]*rho[0])*volFac+((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(403221.4280020346*rdxCp2[0]*rdxCp2Sq[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(93600.0*rdxCp2[0]*rdxCp2Sq[1]+12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-12600.0*rdxCp2[0]*rdxCp2Sq[1])-8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]-75600.0*phiLy[1]*rdxCp2R3[1]+((-82980.0*rdxCp2[0]*phiLy[1])+210000.0*rdxCp2[0]*phiLx[1]-1454922.678357857*rdxCp2[0]*bcVals[1]+(81059.97779422344*phiLy[0]+181865.3347947321*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-11340.0*rdxCp2Sq[0]*phiLy[1])+341400.0*rdxCp2Sq[0]*phiLx[1]-1313587.332460236*rdxCp2Sq[0]*bcVals[1]+(10911.92008768392*phiLy[0]+295661.0728520073*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+50400.0*rdxCp2R3[0]*phiLx[1]-174590.7214029428*rdxCp2R3[0]*bcVals[1]+43647.6803507357*phiLx[0]*rdxCp2R3[0]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+7274.613391789284*rho[0]*rdxCp2Sq[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiLx[3]+((-174590.7214029428*rdxCp2R3[1])-1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]-1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(50400.0*rdxCp2R3[1]+341400.0*rdxCp2[0]*rdxCp2Sq[1]+210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-11340.0*rdxCp2[0]*rdxCp2Sq[1])-82980.0*rdxCp2Sq[0]*rdxCp2[1]-75600.0*rdxCp2R3[0])*phiLx[2]+43647.6803507357*phiLy[0]*rdxCp2R3[1]+((-8640.0*rdxCp2[0]*phiLy[1])+12600.0*rdxCp2[0]*phiLx[1]+87295.36070147139*rdxCp2[0]*bcVals[1]+(295661.0728520073*phiLy[0]+10911.92008768392*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-12600.0*rdxCp2Sq[0]*phiLy[1])+93600.0*rdxCp2Sq[0]*phiLx[1]+403221.4280020346*rdxCp2Sq[0]*bcVals[1]+(181865.3347947321*phiLy[0]+81059.97779422344*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-640.8587988004846*rdxCp2[0]*rdxCp2[1])-2424.871130596428*rdxCp2Sq[0])*rho[2]+((-2424.871130596428*rdxCp2Sq[1])-640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]+5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiLx[3]+(26400.0*rdxCp2[0]*rdxCp2Sq[1]-168000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(16627.68775266122*rdxCp2[0]*rdxCp2Sq[1]+24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-21823.84017536785*rdxCp2[0]*rdxCp2Sq[1])-128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]-29098.45356715714*rdxCp2R3[0])*phiLx[2]-29098.45356715714*phiLy[1]*rdxCp2R3[1]+((-128933.8621154272*rdxCp2[0]*phiLy[1])+24248.71130596428*rdxCp2[0]*phiLx[1]-168000.0*rdxCp2[0]*bcVals[1]+(14400.0*phiLy[0]+21000.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-21823.84017536785*rdxCp2Sq[0]*phiLy[1])+16627.68775266122*rdxCp2Sq[0]*phiLx[1]+26400.0*rdxCp2Sq[0]*bcVals[1]+(21000.0*phiLy[0]+14400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-207.8460969082653*rdxCp2R3[1])-6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]-13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1143.153532995459*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+1091.192008768392*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*volFac+(2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3200.0*rdxCp2R4[1])-96720.0*rdxCp2[0]*rdxCp2R3[1]-222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-1385.640646055102*rdxCp2R4[1])-42816.29596310264*rdxCp2[0]*rdxCp2R3[1]-122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]-72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-623.5382907247956*rdxCp2[0]*rdxCp2R3[1])-22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]-48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-1200.0*phiLy[0]*rdxCp2R4[1]+(2078.460969082652*rdxCp2[0]*phiLy[1]-2078.460969082652*rdxCp2[0]*phiLx[1]-14400.0*rdxCp2[0]*bcVals[1]+((-37080.0*phiLy[0])-1800.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(6131.459858793824*rdxCp2Sq[0]*phiLy[1]-74720.67183852136*rdxCp2Sq[0]*phiLx[1]-359280.0*rdxCp2Sq[0]*bcVals[1]+((-106140.0*phiLy[0])-64710.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(4364.768035073569*rdxCp2R3[0]*phiLy[1]-188308.5637988883*rdxCp2R3[0]*phiLx[1]-879840.0*rdxCp2R3[0]*bcVals[1]+((-63000.0*phiLy[0])-163080.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-65471.52052610354*rdxCp2R4[0]*phiLx[1]-302400.0*rdxCp2R4[0]*bcVals[1]-56700.0*phiLx[0]*rdxCp2R4[0]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = (((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]-11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-10911.92008768392*rdxCp2R3[0]*rho[0])*volFac+(2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-31869.73485926734*rdxCp2[0]*rdxCp2R3[1])-81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]-14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-12000.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25200.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-3600.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+2400.0*phiLy[1]*rdxCp2R4[1]+(24960.0*rdxCp2[0]*phiLy[1]-12000.0*rdxCp2[0]*phiLx[1]+83138.4387633061*rdxCp2[0]*bcVals[1]+((-10392.30484541326*phiLy[0])-10392.30484541326*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(67140.0*rdxCp2Sq[0]*phiLy[1]-114600.0*rdxCp2Sq[0]*phiLx[1]+519615.2422706631*rdxCp2Sq[0]*bcVals[1]+((-30657.29929396912*phiLy[0])-99246.51127369665*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(22680.0*rdxCp2R3[0]*phiLy[1]-237600.0*rdxCp2R3[0]*phiLx[1]+910365.9044582016*rdxCp2R3[0]*bcVals[1]+((-21823.84017536785*phiLy[0])-205767.6359391825*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-75600.0*rdxCp2R4[0]*phiLx[1]+261886.0821044141*rdxCp2R4[0]*bcVals[1]-65471.52052610354*phiLx[0]*rdxCp2R4[0])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+((-120.0*rdxCp2R3[1])-3870.0*rdxCp2[0]*rdxCp2Sq[1]-15150.0*rdxCp2Sq[0]*rdxCp2[1]-15750.0*rdxCp2R3[0])*rho[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+(692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[3]+((-1385.640646055102*rdxCp2R4[1])-43647.6803507357*rdxCp2[0]*rdxCp2R3[1]-145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]-121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(1800.0*rdxCp2[0]*rdxCp2R3[1]+50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+105000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]-50760.0*rdxCp2R3[0]*rdxCp2[1]-56700.0*rdxCp2R4[0])*phiLx[2]+(600.0*rdxCp2[0]*phiLy[1]-600.0*rdxCp2[0]*phiLx[1]-4156.921938165305*rdxCp2[0]*bcVals[1]+(1558.845726811989*phiLy[0]-519.6152422706631*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-21600.0*rdxCp2Sq[0]*phiLx[1])-99766.1265159673*rdxCp2Sq[0]*bcVals[1]+(43647.6803507357*phiLy[0]-18706.14872174387*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-6300.0*rdxCp2R3[0]*phiLy[1])-46800.0*rdxCp2R3[0]*phiLx[1]-201610.7140010173*rdxCp2R3[0]*bcVals[1]+(90932.66739736605*phiLy[0]-40529.98889711172*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiLx[3]+((-10400.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(36373.06695894642*rdxCp2R3[0]*rdxCp2[1]-3464.101615137754*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-2078.460969082652*rdxCp2[0]*rdxCp2R3[1])-39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]-145907.9600296021*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[2]+((-17320.50807568877*rdxCp2[0]*phiLy[1])-3464.101615137754*rdxCp2[0]*phiLx[1]+24000.0*rdxCp2[0]*bcVals[1]+((-3000.0*phiLy[0])-3000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-87295.36070147139*rdxCp2Sq[0]*phiLy[1])-24941.53162899183*rdxCp2Sq[0]*phiLx[1]+86400.0*rdxCp2Sq[0]*bcVals[1]-21600.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-32735.76026305177*rdxCp2R3[0]*phiLy[1])-24941.53162899183*rdxCp2R3[0]*phiLx[1]-39600.0*rdxCp2R3[0]*bcVals[1]+(31500.0*phiLy[0]-21600.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1091.192008768392*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-13146.26562944778*rdxCp2[0]*rdxCp2Sq[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]-207.8460969082653*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*volFac+((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-302400.0*rdxCp2R4[1])-879840.0*rdxCp2[0]*rdxCp2R3[1]-359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-65471.52052610354*rdxCp2R4[1])-188308.5637988883*rdxCp2[0]*rdxCp2R3[1]-74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(4364.768035073569*rdxCp2[0]*rdxCp2R3[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-56700.0*phiLy[0]*rdxCp2R4[1]+((-48635.98667653406*rdxCp2[0]*phiLy[1])-72746.13391789283*rdxCp2[0]*phiLx[1]-42000.0*rdxCp2[0]*bcVals[1]+((-163080.0*phiLy[0])-63000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-22447.37846609264*rdxCp2Sq[0]*phiLy[1])-122559.9151435737*rdxCp2Sq[0]*phiLx[1]-222560.0*rdxCp2Sq[0]*bcVals[1]+((-64710.0*phiLy[0])-106140.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-623.5382907247956*rdxCp2R3[0]*phiLy[1])-42816.29596310264*rdxCp2R3[0]*phiLx[1]-96720.0*rdxCp2R3[0]*bcVals[1]+((-1800.0*phiLy[0])-37080.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-1385.640646055102*rdxCp2R4[0]*phiLx[1]-3200.0*rdxCp2R4[0]*bcVals[1]-1200.0*phiLx[0]*rdxCp2R4[0]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(885.0*rdxCp2[0]*rdxCp2Sq[1]+360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-15750.0*rdxCp2R3[1])-15150.0*rdxCp2[0]*rdxCp2Sq[1]-3870.0*rdxCp2Sq[0]*rdxCp2[1]-120.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*volFac+((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-46800.0*rdxCp2[0]*rdxCp2R3[1])-21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(600.0*rdxCp2R3[0]*rdxCp2[1]-6300.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]-56700.0*phiLy[1]*rdxCp2R4[1]+((-50760.0*rdxCp2[0]*phiLy[1])+105000.0*rdxCp2[0]*phiLx[1]-121243.5565298214*rdxCp2[0]*bcVals[1]+(90932.66739736605*phiLx[0]-40529.98889711172*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-12870.0*rdxCp2Sq[0]*phiLy[1])+50400.0*rdxCp2Sq[0]*phiLx[1]-145838.6779972995*rdxCp2Sq[0]*bcVals[1]+(43647.6803507357*phiLx[0]-18706.14872174387*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-360.0*rdxCp2R3[0]*phiLy[1])+1800.0*rdxCp2R3[0]*phiLx[1]-43647.6803507357*rdxCp2R3[0]*bcVals[1]+(1558.845726811989*phiLx[0]-519.6152422706631*phiLy[0])*rdxCp2R3[0])*rdxCp2[1]-1385.640646055102*rdxCp2R4[0]*bcVals[1]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = (((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-10911.92008768392*rho[0]*rdxCp2R3[1]-30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiLx[3]+(261886.0821044141*rdxCp2R4[1]+910365.9044582016*rdxCp2[0]*rdxCp2R3[1]+519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]+83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-75600.0*rdxCp2R4[1])-237600.0*rdxCp2[0]*rdxCp2R3[1]-114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-12000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(22680.0*rdxCp2[0]*rdxCp2R3[1]+67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiLx[2]-65471.52052610354*phiLy[0]*rdxCp2R4[1]+((-25920.0*rdxCp2[0]*phiLy[1])-25200.0*rdxCp2[0]*phiLx[1]-14549.22678357857*rdxCp2[0]*bcVals[1]+((-205767.6359391825*phiLy[0])-21823.84017536785*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-25920.0*rdxCp2Sq[0]*phiLy[1])-35400.0*rdxCp2Sq[0]*phiLx[1]-81752.798117251*rdxCp2Sq[0]*bcVals[1]+((-99246.51127369665*phiLy[0])-30657.29929396912*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3600.0*rdxCp2R3[0]*phiLy[1])-12000.0*rdxCp2R3[0]*phiLx[1]-31869.73485926734*rdxCp2R3[0]*bcVals[1]+((-10392.30484541326*phiLy[0])-10392.30484541326*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*volFac+((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-39600.0*rdxCp2[0]*rdxCp2R3[1])+86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-24941.53162899183*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-32735.76026305177*rdxCp2[0]*rdxCp2R3[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]-17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-65471.52052610354*phiLy[1]*rdxCp2R4[1]+((-145907.9600296021*rdxCp2[0]*phiLy[1])+36373.06695894642*rdxCp2[0]*phiLx[1]-42000.0*rdxCp2[0]*bcVals[1]+(31500.0*phiLx[0]-21600.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-39386.83536411626*rdxCp2Sq[0]*phiLy[1])-35400.0*rdxCp2Sq[0]*bcVals[1]-21600.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2078.460969082652*rdxCp2R3[0]*phiLy[1])-3464.101615137754*rdxCp2R3[0]*phiLx[1]-10400.0*rdxCp2R3[0]*bcVals[1]+((-3000.0*phiLy[0])-3000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(25.98076211353316*rdxCp2Sq[1]+285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+(285.7883832488647*rdxCp2[0]*rdxCp2[1]+25.98076211353316*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*volFac+(600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(400.0*rdxCp2R3[1]+4440.0*rdxCp2[0]*rdxCp2Sq[1]+200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(173.2050807568877*rdxCp2R3[1]+1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+150.0*phiLy[0]*rdxCp2R3[1]+(519.6152422706631*rdxCp2[0]*phiLy[1]+346.4101615137754*rdxCp2[0]*phiLx[1]+200.0*rdxCp2[0]*bcVals[1]+(1710.0*phiLy[0]+300.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(103.9230484541326*rdxCp2Sq[0]*phiLy[1]+1974.53792062852*rdxCp2Sq[0]*phiLx[1]+4440.0*rdxCp2Sq[0]*bcVals[1]+(300.0*phiLy[0]+1710.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+173.2050807568877*rdxCp2R3[0]*phiLx[1]+400.0*rdxCp2R3[0]*bcVals[1]+150.0*phiLx[0]*rdxCp2R3[0])/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = (((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0]*rho[0])*volFac+(1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(7967.433714816835*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(3000.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(900.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+900.0*phiLy[1]*rdxCp2R3[1]+(3060.0*rdxCp2[0]*phiLy[1]-3000.0*rdxCp2[0]*phiLx[1]+3464.101615137754*rdxCp2[0]*bcVals[1]+(2598.076211353316*phiLy[0]-2598.076211353316*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(360.0*rdxCp2Sq[0]*phiLy[1]-600.0*rdxCp2Sq[0]*phiLx[1]+12124.35565298214*rdxCp2Sq[0]*bcVals[1]+(519.6152422706631*phiLy[0]-519.6152422706631*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+1039.230484541326*rdxCp2R3[0]*bcVals[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = (((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+259.8076211353315*rho[0]*rdxCp2Sq[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+(1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiLx[3]+(1039.230484541326*rdxCp2R3[1]+12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+3060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiLx[2]+(900.0*rdxCp2[0]*phiLy[1]+600.0*rdxCp2[0]*phiLx[1]+346.4101615137754*rdxCp2[0]*bcVals[1]+(519.6152422706631*phiLx[0]-519.6152422706631*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-900.0*rdxCp2Sq[0]*phiLy[1])+3000.0*rdxCp2Sq[0]*phiLx[1]+7967.433714816835*rdxCp2Sq[0]*bcVals[1]+(2598.076211353316*phiLx[0]-2598.076211353316*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0])*rho[2]+(129.9038105676658*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*volFac+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(2600.0*rdxCp2[0]*rdxCp2Sq[1]+1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(866.0254037844386*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2Sq[1])-2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-2598.076211353316*rdxCp2[0]*phiLy[1])-866.0254037844386*rdxCp2[0]*phiLx[1]+1000.0*rdxCp2[0]*bcVals[1]+(750.0*phiLy[0]-750.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-519.6152422706631*rdxCp2Sq[0]*phiLy[1])+866.0254037844386*rdxCp2Sq[0]*phiLx[1]+2600.0*rdxCp2Sq[0]*bcVals[1]+(750.0*phiLx[0]-750.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = 1.0/dxLx[1]; 
  rdxUx[1]   = 1.0/dxUx[1]; 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = 1.0/dxLy[0]; 
  rdxUy[0]   = 1.0/dxUy[0]; 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = 1.0/dxLy[1]; 
  rdxUy[1]   = 1.0/dxUy[1]; 
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

  phiC[0] = ((((24600.0*rdxUx[0]-24600.0*rdxLx[0])*rdxUySq[1]+(24600.0*rdxUxSq[0]-24600.0*rdxLxSq[0])*rdxUy[1]+(24600.0*rdxLx[0]-24600.0*rdxUx[0])*rdxLySq[1]+(24600.0*rdxLxSq[0]-24600.0*rdxUxSq[0])*rdxLy[1])*rho[3]+((-1662.768775266122*rdxUyCu[1])+((-109742.739167564*rdxLy[1])-65332.95646149805*rdxUx[0]-65332.95646149805*rdxLx[0])*rdxUySq[1]+(109742.739167564*rdxLySq[1]-63670.18768623193*rdxUxSq[0]-19260.40498016591*rdxLx[0]*rdxUx[0]-63670.18768623193*rdxLxSq[0])*rdxUy[1]+1662.768775266122*rdxLyCu[1]+(65332.95646149805*rdxUx[0]+65332.95646149805*rdxLx[0])*rdxLySq[1]+(63670.18768623193*rdxUxSq[0]+19260.40498016591*rdxLx[0]*rdxUx[0]+63670.18768623193*rdxLxSq[0])*rdxLy[1])*rho[2]+((63670.18768623193*rdxLx[0]-63670.18768623193*rdxUx[0])*rdxUySq[1]+((19260.40498016591*rdxLx[0]-19260.40498016591*rdxUx[0])*rdxLy[1]-65332.95646149805*rdxUxSq[0]+65332.95646149805*rdxLxSq[0])*rdxUy[1]+(63670.18768623193*rdxLx[0]-63670.18768623193*rdxUx[0])*rdxLySq[1]+(65332.95646149805*rdxLxSq[0]-65332.95646149805*rdxUxSq[0])*rdxLy[1]-1662.768775266122*rdxUxCu[0]-109742.739167564*rdxLx[0]*rdxUxSq[0]+109742.739167564*rdxLxSq[0]*rdxUx[0]+1662.768775266122*rdxLxCu[0])*rho[1]+4416.0*rho[0]*rdxUyCu[1]+(300288.0*rho[0]*rdxLy[1]+(176968.0*rdxUx[0]+176968.0*rdxLx[0])*rho[0])*rdxUySq[1]+(300288.0*rho[0]*rdxLySq[1]+(578576.0*rdxUx[0]+578576.0*rdxLx[0])*rho[0]*rdxLy[1]+(176968.0*rdxUxSq[0]+578576.0*rdxLx[0]*rdxUx[0]+176968.0*rdxLxSq[0])*rho[0])*rdxUy[1]+4416.0*rho[0]*rdxLyCu[1]+(176968.0*rdxUx[0]+176968.0*rdxLx[0])*rho[0]*rdxLySq[1]+(176968.0*rdxUxSq[0]+578576.0*rdxLx[0]*rdxUx[0]+176968.0*rdxLxSq[0])*rho[0]*rdxLy[1]+(4416.0*rdxUxCu[0]+300288.0*rdxLx[0]*rdxUxSq[0]+300288.0*rdxLxSq[0]*rdxUx[0]+4416.0*rdxLxCu[0])*rho[0])*omega*volFac+(((47400.0*rdxUx[0]-47400.0*rdxLx[0])*rdxUyCu[1]+((20850.0*rdxUx[0]-20850.0*rdxLx[0])*rdxLy[1]+49200.0*rdxUxSq[0]-49200.0*rdxLxSq[0])*rdxUySq[1]+((90450.0*rdxUx[0]-90450.0*rdxLx[0])*rdxLySq[1]+(92250.0*rdxUxSq[0]-92250.0*rdxLxSq[0])*rdxLy[1]+1800.0*rdxUxCu[0]+118800.0*rdxLx[0]*rdxUxSq[0]-118800.0*rdxLxSq[0]*rdxUx[0]-1800.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(1800.0*rdxUx[0]*rdxUyCu[1]+(118800.0*rdxUx[0]*rdxLy[1]+49200.0*rdxUxSq[0]+92250.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-118800.0*rdxUx[0]*rdxLySq[1])+47400.0*rdxUxCu[0]+20850.0*rdxLx[0]*rdxUxSq[0]+90450.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-1800.0*rdxUx[0]*rdxLyCu[1]+((-49200.0*rdxUxSq[0])-92250.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-47400.0*rdxUxCu[0])-20850.0*rdxLx[0]*rdxUxSq[0]-90450.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((90450.0*rdxLx[0]-90450.0*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((20850.0*rdxLx[0]-20850.0*rdxUx[0])*rdxLySq[1]+(92250.0*rdxLxSq[0]-92250.0*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(47400.0*rdxLx[0]-47400.0*rdxUx[0])*rdxLyCu[1]+(49200.0*rdxLxSq[0]-49200.0*rdxUxSq[0])*rdxLySq[1]+((-1800.0*rdxUxCu[0])-118800.0*rdxLx[0]*rdxUxSq[0]+118800.0*rdxLxSq[0]*rdxUx[0]+1800.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-1800.0*rdxLx[0]*rdxUyCu[1])+((-118800.0*rdxLx[0]*rdxLy[1])-92250.0*rdxLx[0]*rdxUx[0]-49200.0*rdxLxSq[0])*rdxUySq[1]+(118800.0*rdxLx[0]*rdxLySq[1]-90450.0*rdxLx[0]*rdxUxSq[0]-20850.0*rdxLxSq[0]*rdxUx[0]-47400.0*rdxLxCu[0])*rdxUy[1]+1800.0*rdxLx[0]*rdxLyCu[1]+(92250.0*rdxLx[0]*rdxUx[0]+49200.0*rdxLxSq[0])*rdxLySq[1]+(90450.0*rdxLx[0]*rdxUxSq[0]+20850.0*rdxLxSq[0]*rdxUx[0]+47400.0*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((-1662.768775266122*rdxUyR4[1])+((-114523.1993964541*rdxLy[1])-67203.57133367243*rdxUx[0]-67203.57133367243*rdxLx[0])*rdxUyCu[1]+((-210548.0961680727*rdxLySq[1])+((-313163.4462624908*rdxUx[0])-313163.4462624908*rdxLx[0])*rdxLy[1]-67931.03267285136*rdxUxSq[0]-304737.0190836682*rdxLx[0]*rdxUx[0]-67931.03267285136*rdxLxSq[0])*rdxUySq[1]+((-3117.691453623978*rdxLyCu[1])+((-124369.9082374832*rdxUx[0])-124369.9082374832*rdxLx[0])*rdxLySq[1]+((-123642.4468983043*rdxUxSq[0])-321589.8734413133*rdxLx[0]*rdxUx[0]-123642.4468983043*rdxLxSq[0])*rdxLy[1]-2390.23011444505*rdxUxCu[0]-162535.6477822634*rdxLx[0]*rdxUxSq[0]-162535.6477822634*rdxLxSq[0]*rdxUx[0]-2390.23011444505*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-1870.614872174387*rdxUx[0]*rdxUyCu[1])+((-123460.5815635095*rdxUx[0]*rdxLy[1])-46869.29485281381*rdxUxSq[0]-100129.8571855568*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(123460.5815635095*rdxUx[0]*rdxLySq[1]-44998.67998063943*rdxUxCu[0]-21667.95560268665*rdxLx[0]*rdxUxSq[0]-98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+1870.614872174387*rdxUx[0]*rdxLyCu[1]+(46869.29485281381*rdxUxSq[0]+100129.8571855568*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(44998.67998063943*rdxUxCu[0]+21667.95560268665*rdxLx[0]*rdxUxSq[0]+98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+(3117.691453623978*rdxLy[1]*rdxUyCu[1]+(210548.0961680727*rdxLySq[1]+(124369.9082374832*rdxUx[0]+124369.9082374832*rdxLx[0])*rdxLy[1])*rdxUySq[1]+(114523.1993964541*rdxLyCu[1]+(313163.4462624908*rdxUx[0]+313163.4462624908*rdxLx[0])*rdxLySq[1]+(123642.4468983043*rdxUxSq[0]+321589.8734413133*rdxLx[0]*rdxUx[0]+123642.4468983043*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+1662.768775266122*rdxLyR4[1]+(67203.57133367243*rdxUx[0]+67203.57133367243*rdxLx[0])*rdxLyCu[1]+(67931.03267285136*rdxUxSq[0]+304737.0190836682*rdxLx[0]*rdxUx[0]+67931.03267285136*rdxLxSq[0])*rdxLySq[1]+(2390.23011444505*rdxUxCu[0]+162535.6477822634*rdxLx[0]*rdxUxSq[0]+162535.6477822634*rdxLxSq[0]*rdxUx[0]+2390.23011444505*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-1870.614872174387*rdxLx[0]*rdxUyCu[1])+((-123460.5815635095*rdxLx[0]*rdxLy[1])-100129.8571855568*rdxLx[0]*rdxUx[0]-46869.29485281381*rdxLxSq[0])*rdxUySq[1]+(123460.5815635095*rdxLx[0]*rdxLySq[1]-98259.24231338239*rdxLx[0]*rdxUxSq[0]-21667.95560268665*rdxLxSq[0]*rdxUx[0]-44998.67998063943*rdxLxCu[0])*rdxUy[1]+1870.614872174387*rdxLx[0]*rdxLyCu[1]+(100129.8571855568*rdxLx[0]*rdxUx[0]+46869.29485281381*rdxLxSq[0])*rdxLySq[1]+(98259.24231338239*rdxLx[0]*rdxUxSq[0]+21667.95560268665*rdxLxSq[0]*rdxUx[0]+44998.67998063943*rdxLxCu[0])*rdxLy[1])*phiLx[2]+(1584.0*phiUy[0]-144.0*phiC[0])*rdxUyR4[1]+((109512.0*phiUy[0]+3384.0*phiLy[0]-19296.0*phiC[0])*rdxLy[1]+(44998.67998063943*rdxLx[0]-44998.67998063943*rdxUx[0])*phiUy[1]-2390.23011444505*rdxUx[0]*phiUx[1]+2390.23011444505*rdxLx[0]*phiLx[1]+(64182.0*phiUy[0]+2484.0*phiUx[0]-10086.0*phiC[0])*rdxUx[0]+(64182.0*phiUy[0]+2484.0*phiLx[0]-10086.0*phiC[0])*rdxLx[0])*rdxUyCu[1]+((228312.0*phiUy[0]+228312.0*phiLy[0]-646704.0*phiC[0])*rdxLySq[1]+((21667.95560268665*rdxLx[0]-21667.95560268665*rdxUx[0])*phiUy[1]-162535.6477822634*rdxUx[0]*phiUx[1]+(98259.24231338239*rdxLx[0]-98259.24231338239*rdxUx[0])*phiLy[1]+162535.6477822634*rdxLx[0]*phiLx[1]+(325449.0*phiUy[0]+168912.0*phiUx[0]+134907.0*phiLy[0]-685848.0*phiC[0])*rdxUx[0]+(325449.0*phiUy[0]+134907.0*phiLy[0]+168912.0*phiLx[0]-685848.0*phiC[0])*rdxLx[0])*rdxLy[1]+(46869.29485281381*rdxLxSq[0]-46869.29485281381*rdxUxSq[0])*phiUy[1]+((-67931.03267285136*rdxUxSq[0])-123642.4468983043*rdxLx[0]*rdxUx[0])*phiUx[1]+(123642.4468983043*rdxLx[0]*rdxUx[0]+67931.03267285136*rdxLxSq[0])*phiLx[1]+(65082.0*phiUy[0]+65082.0*phiUx[0]-19884.0*phiC[0])*rdxUxSq[0]+(315024.0*phiUy[0]+134007.0*phiUx[0]+134007.0*phiLx[0]-676638.0*phiC[0])*rdxLx[0]*rdxUx[0]+(65082.0*phiUy[0]+65082.0*phiLx[0]-19884.0*phiC[0])*rdxLxSq[0])*rdxUySq[1]+((3384.0*phiUy[0]+109512.0*phiLy[0]-19296.0*phiC[0])*rdxLyCu[1]+((98259.24231338239*rdxLx[0]-98259.24231338239*rdxUx[0])*phiUy[1]-162535.6477822634*rdxUx[0]*phiUx[1]+(21667.95560268665*rdxLx[0]-21667.95560268665*rdxUx[0])*phiLy[1]+162535.6477822634*rdxLx[0]*phiLx[1]+(134907.0*phiUy[0]+168912.0*phiUx[0]+325449.0*phiLy[0]-685848.0*phiC[0])*rdxUx[0]+(134907.0*phiUy[0]+325449.0*phiLy[0]+168912.0*phiLx[0]-685848.0*phiC[0])*rdxLx[0])*rdxLySq[1]+((100129.8571855568*rdxLxSq[0]-100129.8571855568*rdxUxSq[0])*phiUy[1]+((-304737.0190836682*rdxUxSq[0])-321589.8734413133*rdxLx[0]*rdxUx[0])*phiUx[1]+(100129.8571855568*rdxLxSq[0]-100129.8571855568*rdxUxSq[0])*phiLy[1]+(321589.8734413133*rdxLx[0]*rdxUx[0]+304737.0190836682*rdxLxSq[0])*phiLx[1]+(134007.0*phiUy[0]+315024.0*phiUx[0]+134007.0*phiLy[0]-676638.0*phiC[0])*rdxUxSq[0]+(335874.0*phiUy[0]+335874.0*phiUx[0]+335874.0*phiLy[0]+335874.0*phiLx[0]-1410216.0*phiC[0])*rdxLx[0]*rdxUx[0]+(134007.0*phiUy[0]+134007.0*phiLy[0]+315024.0*phiLx[0]-676638.0*phiC[0])*rdxLxSq[0])*rdxLy[1]+((-1870.614872174387*rdxUxCu[0])-123460.5815635095*rdxLx[0]*rdxUxSq[0]+123460.5815635095*rdxLxSq[0]*rdxUx[0]+1870.614872174387*rdxLxCu[0])*phiUy[1]+((-67203.57133367243*rdxUxCu[0])-313163.4462624908*rdxLx[0]*rdxUxSq[0]-124369.9082374832*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(124369.9082374832*rdxLx[0]*rdxUxSq[0]+313163.4462624908*rdxLxSq[0]*rdxUx[0]+67203.57133367243*rdxLxCu[0])*phiLx[1]+(2484.0*phiUy[0]+64182.0*phiUx[0]-10086.0*phiC[0])*rdxUxCu[0]+(168912.0*phiUy[0]+325449.0*phiUx[0]+134907.0*phiLx[0]-685848.0*phiC[0])*rdxLx[0]*rdxUxSq[0]+(168912.0*phiUy[0]+134907.0*phiUx[0]+325449.0*phiLx[0]-685848.0*phiC[0])*rdxLxSq[0]*rdxUx[0]+(2484.0*phiUy[0]+64182.0*phiLx[0]-10086.0*phiC[0])*rdxLxCu[0])*rdxUy[1]+(1584.0*phiLy[0]-144.0*phiC[0])*rdxLyR4[1]+((-2390.23011444505*rdxUx[0]*phiUx[1])+(44998.67998063943*rdxLx[0]-44998.67998063943*rdxUx[0])*phiLy[1]+2390.23011444505*rdxLx[0]*phiLx[1]+(2484.0*phiUx[0]+64182.0*phiLy[0]-10086.0*phiC[0])*rdxUx[0]+(64182.0*phiLy[0]+2484.0*phiLx[0]-10086.0*phiC[0])*rdxLx[0])*rdxLyCu[1]+(((-67931.03267285136*rdxUxSq[0])-123642.4468983043*rdxLx[0]*rdxUx[0])*phiUx[1]+(46869.29485281381*rdxLxSq[0]-46869.29485281381*rdxUxSq[0])*phiLy[1]+(123642.4468983043*rdxLx[0]*rdxUx[0]+67931.03267285136*rdxLxSq[0])*phiLx[1]+(65082.0*phiUx[0]+65082.0*phiLy[0]-19884.0*phiC[0])*rdxUxSq[0]+(134007.0*phiUx[0]+315024.0*phiLy[0]+134007.0*phiLx[0]-676638.0*phiC[0])*rdxLx[0]*rdxUx[0]+(65082.0*phiLy[0]+65082.0*phiLx[0]-19884.0*phiC[0])*rdxLxSq[0])*rdxLySq[1]+(((-67203.57133367243*rdxUxCu[0])-313163.4462624908*rdxLx[0]*rdxUxSq[0]-124369.9082374832*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-1870.614872174387*rdxUxCu[0])-123460.5815635095*rdxLx[0]*rdxUxSq[0]+123460.5815635095*rdxLxSq[0]*rdxUx[0]+1870.614872174387*rdxLxCu[0])*phiLy[1]+(124369.9082374832*rdxLx[0]*rdxUxSq[0]+313163.4462624908*rdxLxSq[0]*rdxUx[0]+67203.57133367243*rdxLxCu[0])*phiLx[1]+(64182.0*phiUx[0]+2484.0*phiLy[0]-10086.0*phiC[0])*rdxUxCu[0]+(325449.0*phiUx[0]+168912.0*phiLy[0]+134907.0*phiLx[0]-685848.0*phiC[0])*rdxLx[0]*rdxUxSq[0]+(134907.0*phiUx[0]+168912.0*phiLy[0]+325449.0*phiLx[0]-685848.0*phiC[0])*rdxLxSq[0]*rdxUx[0]+(2484.0*phiLy[0]+64182.0*phiLx[0]-10086.0*phiC[0])*rdxLxCu[0])*rdxLy[1]+((-1662.768775266122*rdxUxR4[0])-114523.1993964541*rdxLx[0]*rdxUxCu[0]-210548.0961680727*rdxLxSq[0]*rdxUxSq[0]-3117.691453623978*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(3117.691453623978*rdxLx[0]*rdxUxCu[0]+210548.0961680727*rdxLxSq[0]*rdxUxSq[0]+114523.1993964541*rdxLxCu[0]*rdxUx[0]+1662.768775266122*rdxLxR4[0])*phiLx[1]+(1584.0*phiUx[0]-144.0*phiC[0])*rdxUxR4[0]+(109512.0*phiUx[0]+3384.0*phiLx[0]-19296.0*phiC[0])*rdxLx[0]*rdxUxCu[0]+(228312.0*phiUx[0]+228312.0*phiLx[0]-646704.0*phiC[0])*rdxLxSq[0]*rdxUxSq[0]+(3384.0*phiUx[0]+109512.0*phiLx[0]-19296.0*phiC[0])*rdxLxCu[0]*rdxUx[0]+(1584.0*phiLx[0]-144.0*phiC[0])*rdxLxR4[0])*omega+144.0*phiC[0]*rdxUyR4[1]+(19296.0*phiC[0]*rdxLy[1]+10086.0*phiC[0]*rdxUx[0]+10086.0*phiC[0]*rdxLx[0])*rdxUyCu[1]+(646704.0*phiC[0]*rdxLySq[1]+(685848.0*phiC[0]*rdxUx[0]+685848.0*phiC[0]*rdxLx[0])*rdxLy[1]+19884.0*phiC[0]*rdxUxSq[0]+676638.0*phiC[0]*rdxLx[0]*rdxUx[0]+19884.0*phiC[0]*rdxLxSq[0])*rdxUySq[1]+(19296.0*phiC[0]*rdxLyCu[1]+(685848.0*phiC[0]*rdxUx[0]+685848.0*phiC[0]*rdxLx[0])*rdxLySq[1]+(676638.0*phiC[0]*rdxUxSq[0]+1410216.0*phiC[0]*rdxLx[0]*rdxUx[0]+676638.0*phiC[0]*rdxLxSq[0])*rdxLy[1]+10086.0*phiC[0]*rdxUxCu[0]+685848.0*phiC[0]*rdxLx[0]*rdxUxSq[0]+685848.0*phiC[0]*rdxLxSq[0]*rdxUx[0]+10086.0*phiC[0]*rdxLxCu[0])*rdxUy[1]+144.0*phiC[0]*rdxLyR4[1]+(10086.0*phiC[0]*rdxUx[0]+10086.0*phiC[0]*rdxLx[0])*rdxLyCu[1]+(19884.0*phiC[0]*rdxUxSq[0]+676638.0*phiC[0]*rdxLx[0]*rdxUx[0]+19884.0*phiC[0]*rdxLxSq[0])*rdxLySq[1]+(10086.0*phiC[0]*rdxUxCu[0]+685848.0*phiC[0]*rdxLx[0]*rdxUxSq[0]+685848.0*phiC[0]*rdxLxSq[0]*rdxUx[0]+10086.0*phiC[0]*rdxLxCu[0])*rdxLy[1]+144.0*phiC[0]*rdxUxR4[0]+19296.0*phiC[0]*rdxLx[0]*rdxUxCu[0]+646704.0*phiC[0]*rdxLxSq[0]*rdxUxSq[0]+19296.0*phiC[0]*rdxLxCu[0]*rdxUx[0]+144.0*phiC[0]*rdxLxR4[0])/(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0]); 
  phiC[1] = -(1.0*(((831.384387633061*rdxUyCu[1]+(54871.36958378201*rdxLy[1]+25565.06991971662*rdxUx[0]+25565.06991971662*rdxLx[0])*rdxUySq[1]+((-54871.36958378201*rdxLySq[1])+24733.68553208356*rdxUxSq[0]-4572.614131981835*rdxLx[0]*rdxUx[0]+24733.68553208356*rdxLxSq[0])*rdxUy[1]-831.384387633061*rdxLyCu[1]+((-25565.06991971662*rdxUx[0])-25565.06991971662*rdxLx[0])*rdxLySq[1]+((-24733.68553208356*rdxUxSq[0])+4572.614131981835*rdxLx[0]*rdxUx[0]-24733.68553208356*rdxLxSq[0])*rdxLy[1])*rho[3]+((63960.0*rdxLx[0]-63960.0*rdxUx[0])*rdxUySq[1]+(63960.0*rdxLxSq[0]-63960.0*rdxUxSq[0])*rdxUy[1]+(63960.0*rdxUx[0]-63960.0*rdxLx[0])*rdxLySq[1]+(63960.0*rdxUxSq[0]-63960.0*rdxLxSq[0])*rdxLy[1])*rho[2]+((-2208.0*rdxUyCu[1])+((-150144.0*rdxLy[1])-70104.0*rdxUx[0]-70104.0*rdxLx[0])*rdxUySq[1]+((-150144.0*rdxLySq[1])+((-283728.0*rdxUx[0])-283728.0*rdxLx[0])*rdxLy[1]-69624.0*rdxUxSq[0]-251568.0*rdxLx[0]*rdxUx[0]-69624.0*rdxLxSq[0])*rdxUy[1]-2208.0*rdxLyCu[1]+((-70104.0*rdxUx[0])-70104.0*rdxLx[0])*rdxLySq[1]+((-69624.0*rdxUxSq[0])-251568.0*rdxLx[0]*rdxUx[0]-69624.0*rdxLxSq[0])*rdxLy[1]-1728.0*rdxUxCu[0]-117504.0*rdxLx[0]*rdxUxSq[0]-117504.0*rdxLxSq[0]*rdxUx[0]-1728.0*rdxLxCu[0])*rho[1]+(165542.487984203*rdxUx[0]-165542.487984203*rdxLx[0])*rho[0]*rdxUySq[1]+((50077.05294843137*rdxUx[0]-50077.05294843137*rdxLx[0])*rho[0]*rdxLy[1]+(169865.6867998949*rdxUxSq[0]-169865.6867998949*rdxLxSq[0])*rho[0])*rdxUy[1]+(165542.487984203*rdxUx[0]-165542.487984203*rdxLx[0])*rho[0]*rdxLySq[1]+(169865.6867998949*rdxUxSq[0]-169865.6867998949*rdxLxSq[0])*rho[0]*rdxLy[1]+(4323.198815691917*rdxUxCu[0]+285331.1218356665*rdxLx[0]*rdxUxSq[0]-285331.1218356665*rdxLxSq[0]*rdxUx[0]-4323.198815691917*rdxLxCu[0])*rho[0])*omega*volFac+((1662.768775266122*rdxUyR4[1]+(114523.1993964541*rdxLy[1]+53520.3699538783*rdxUx[0]+53520.3699538783*rdxLx[0])*rdxUyCu[1]+(210548.0961680727*rdxLySq[1]+(307144.569706189*rdxUx[0]+307144.569706189*rdxLx[0])*rdxLy[1]+53728.21605078656*rdxUxSq[0]+276331.3858395386*rdxLx[0]*rdxUx[0]+53728.21605078656*rdxLxSq[0])*rdxUySq[1]+(3117.691453623978*rdxLyCu[1]+(98259.24231338239*rdxUx[0]+98259.24231338239*rdxLx[0])*rdxLySq[1]+(97012.16573193281*rdxUxSq[0]+268329.3111085704*rdxLx[0]*rdxUx[0]+97012.16573193281*rdxLxSq[0])*rdxLy[1]+1870.614872174387*rdxUxCu[0]+127201.8113078583*rdxLx[0]*rdxUxSq[0]+127201.8113078583*rdxLxSq[0]*rdxUx[0]+1870.614872174387*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-727.4613391789284*rdxUx[0]*rdxUyCu[1])+((-48012.44838580926*rdxUx[0]*rdxLy[1])+46869.29485281381*rdxUxSq[0]-91608.1672123179*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(48012.44838580926*rdxUx[0]*rdxLySq[1]+47596.75619199274*rdxUxCu[0]+4001.037365484106*rdxLx[0]*rdxUxSq[0]-90880.70587313897*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+727.4613391789284*rdxUx[0]*rdxLyCu[1]+(91608.1672123179*rdxLx[0]*rdxUx[0]-46869.29485281381*rdxUxSq[0])*rdxLySq[1]+((-47596.75619199274*rdxUxCu[0])-4001.037365484106*rdxLx[0]*rdxUxSq[0]+90880.70587313897*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[3]+((-3117.691453623978*rdxLy[1]*rdxUyCu[1])+(((-98259.24231338239*rdxUx[0])-98259.24231338239*rdxLx[0])*rdxLy[1]-210548.0961680727*rdxLySq[1])*rdxUySq[1]+((-114523.1993964541*rdxLyCu[1])+((-307144.569706189*rdxUx[0])-307144.569706189*rdxLx[0])*rdxLySq[1]+((-97012.16573193281*rdxUxSq[0])-268329.3111085704*rdxLx[0]*rdxUx[0]-97012.16573193281*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-1662.768775266122*rdxLyR4[1]+((-53520.3699538783*rdxUx[0])-53520.3699538783*rdxLx[0])*rdxLyCu[1]+((-53728.21605078656*rdxUxSq[0])-276331.3858395386*rdxLx[0]*rdxUx[0]-53728.21605078656*rdxLxSq[0])*rdxLySq[1]+((-1870.614872174387*rdxUxCu[0])-127201.8113078583*rdxLx[0]*rdxUxSq[0]-127201.8113078583*rdxLxSq[0]*rdxUx[0]-1870.614872174387*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-727.4613391789284*rdxLx[0]*rdxUyCu[1])+((-48012.44838580926*rdxLx[0]*rdxLy[1])-91608.1672123179*rdxLx[0]*rdxUx[0]+46869.29485281381*rdxLxSq[0])*rdxUySq[1]+(48012.44838580926*rdxLx[0]*rdxLySq[1]-90880.70587313897*rdxLx[0]*rdxUxSq[0]+4001.037365484106*rdxLxSq[0]*rdxUx[0]+47596.75619199274*rdxLxCu[0])*rdxUy[1]+727.4613391789284*rdxLx[0]*rdxLyCu[1]+(91608.1672123179*rdxLx[0]*rdxUx[0]-46869.29485281381*rdxLxSq[0])*rdxLySq[1]+(90880.70587313897*rdxLx[0]*rdxUxSq[0]-4001.037365484106*rdxLxSq[0]*rdxUx[0]-47596.75619199274*rdxLxCu[0])*rdxLy[1])*phiLx[3]+((61620.0*rdxLx[0]-61620.0*rdxUx[0])*rdxUyCu[1]+((27105.0*rdxLx[0]-27105.0*rdxUx[0])*rdxLy[1]-63960.0*rdxUxSq[0]+63960.0*rdxLxSq[0])*rdxUySq[1]+((117585.0*rdxLx[0]-117585.0*rdxUx[0])*rdxLySq[1]+(119925.0*rdxLxSq[0]-119925.0*rdxUxSq[0])*rdxLy[1]-2340.0*rdxUxCu[0]-154440.0*rdxLx[0]*rdxUxSq[0]+154440.0*rdxLxSq[0]*rdxUx[0]+2340.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(900.0*rdxUx[0]*rdxUyCu[1]+(59400.0*rdxUx[0]*rdxLy[1]-44280.0*rdxUxSq[0]+99630.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-59400.0*rdxUx[0]*rdxLySq[1])-45180.0*rdxUxCu[0]-4950.0*rdxLx[0]*rdxUxSq[0]+98730.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-900.0*rdxUx[0]*rdxLyCu[1]+(44280.0*rdxUxSq[0]-99630.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(45180.0*rdxUxCu[0]+4950.0*rdxLx[0]*rdxUxSq[0]-98730.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1])*phiUx[2]+((117585.0*rdxUx[0]-117585.0*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((27105.0*rdxUx[0]-27105.0*rdxLx[0])*rdxLySq[1]+(119925.0*rdxUxSq[0]-119925.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(61620.0*rdxUx[0]-61620.0*rdxLx[0])*rdxLyCu[1]+(63960.0*rdxUxSq[0]-63960.0*rdxLxSq[0])*rdxLySq[1]+(2340.0*rdxUxCu[0]+154440.0*rdxLx[0]*rdxUxSq[0]-154440.0*rdxLxSq[0]*rdxUx[0]-2340.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-900.0*rdxLx[0]*rdxUyCu[1])+((-59400.0*rdxLx[0]*rdxLy[1])-99630.0*rdxLx[0]*rdxUx[0]+44280.0*rdxLxSq[0])*rdxUySq[1]+(59400.0*rdxLx[0]*rdxLySq[1]-98730.0*rdxLx[0]*rdxUxSq[0]+4950.0*rdxLxSq[0]*rdxUx[0]+45180.0*rdxLxCu[0])*rdxUy[1]+900.0*rdxLx[0]*rdxLyCu[1]+(99630.0*rdxLx[0]*rdxUx[0]-44280.0*rdxLxSq[0])*rdxLySq[1]+(98730.0*rdxLx[0]*rdxUxSq[0]-4950.0*rdxLxSq[0]*rdxUx[0]-45180.0*rdxLxCu[0])*rdxLy[1])*phiLx[2]+(144.0*phiC[1]-1584.0*phiUy[1])*rdxUyR4[1]+(((-109512.0*phiUy[1])-3384.0*phiLy[1]+19296.0*phiC[1])*rdxLy[1]+((-51192.0*rdxUx[0])-51192.0*rdxLx[0])*phiUy[1]+966.0*rdxUx[0]*phiUx[1]+966.0*rdxLx[0]*phiLx[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*phiC[1]+(58498.28397483125*phiUy[0]-1195.115057222525*phiUx[0])*rdxUx[0]+(1195.115057222525*phiLx[0]-58498.28397483125*phiUy[0])*rdxLx[0])*rdxUyCu[1]+(((-228312.0*phiUy[1])-228312.0*phiLy[1]+646704.0*phiC[1])*rdxLySq[1]+(((-319194.0*rdxUx[0])-319194.0*rdxLx[0])*phiUy[1]+65688.0*rdxUx[0]*phiUx[1]+((-106542.0*rdxUx[0])-106542.0*rdxLx[0])*phiLy[1]+65688.0*rdxLx[0]*phiLx[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*phiC[1]+(28168.34228349264*phiUy[0]-81267.82389113172*phiUx[0]+127737.0150073971*phiLy[0])*rdxUx[0]+((-28168.34228349264*phiUy[0])-127737.0150073971*phiLy[0]+81267.82389113172*phiLx[0])*rdxLx[0])*rdxLy[1]+((-51552.0*rdxUxSq[0])-287964.0*rdxLx[0]*rdxUx[0]-51552.0*rdxLxSq[0])*phiUy[1]+(120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxUxSq[0])*phiUx[1]+(120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*phiLx[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*phiC[1]+(60930.08330865796*phiUy[0]+55172.74642429901*phiUx[0])*rdxUxSq[0]+(131062.5525579294*phiLx[0]-131062.5525579294*phiUx[0])*rdxLx[0]*rdxUx[0]+((-60930.08330865796*phiUy[0])-55172.74642429901*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-3384.0*phiUy[1])-109512.0*phiLy[1]+19296.0*phiC[1])*rdxLyCu[1]+(((-106542.0*rdxUx[0])-106542.0*rdxLx[0])*phiUy[1]+65688.0*rdxUx[0]*phiUx[1]+((-319194.0*rdxUx[0])-319194.0*rdxLx[0])*phiLy[1]+65688.0*rdxLx[0]*phiLx[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*phiC[1]+(127737.0150073971*phiUy[0]-81267.82389113172*phiUx[0]+28168.34228349264*phiLy[0])*rdxUx[0]+((-127737.0150073971*phiUy[0])-28168.34228349264*phiLy[0]+81267.82389113172*phiLx[0])*rdxLx[0])*rdxLySq[1]+(((-105102.0*rdxUxSq[0])-278064.0*rdxLx[0]*rdxUx[0]-105102.0*rdxLxSq[0])*phiUy[1]+(97026.0*rdxUxSq[0]+151236.0*rdxLx[0]*rdxUx[0])*phiUx[1]+((-105102.0*rdxUxSq[0])-278064.0*rdxLx[0]*rdxUx[0]-105102.0*rdxLxSq[0])*phiLy[1]+(151236.0*rdxLx[0]*rdxUx[0]+97026.0*rdxLxSq[0])*phiLx[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*phiC[1]+(130168.8143412238*phiUy[0]-125403.9425696018*phiUx[0]+130168.8143412238*phiLy[0])*rdxUxSq[0]+(181740.6271365871*phiLx[0]-181740.6271365871*phiUx[0])*rdxLx[0]*rdxUx[0]+((-130168.8143412238*phiUy[0])-130168.8143412238*phiLy[0]+125403.9425696018*phiLx[0])*rdxLxSq[0])*rdxLy[1]+((-1944.0*rdxUxCu[0])-132192.0*rdxLx[0]*rdxUxSq[0]-132192.0*rdxLxSq[0]*rdxUx[0]-1944.0*rdxLxCu[0])*phiUy[1]+((-61482.0*rdxUxCu[0])+110061.0*rdxLx[0]*rdxUxSq[0]+122403.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(122403.0*rdxLx[0]*rdxUxSq[0]+110061.0*rdxLxSq[0]*rdxUx[0]-61482.0*rdxLxCu[0])*phiLx[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*phiC[1]+(2431.799333826703*phiUy[0]+57864.35337926103*phiUx[0])*rdxUxCu[0]+(160498.7560325623*phiUy[0]-136165.1742370272*phiUx[0]+133234.5442706207*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-160498.7560325623*phiUy[0])-133234.5442706207*phiUx[0]+136165.1742370272*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-2431.799333826703*phiUy[0])-57864.35337926103*phiLx[0])*rdxLxCu[0])*rdxUy[1]+(144.0*phiC[1]-1584.0*phiLy[1])*rdxLyR4[1]+(966.0*rdxUx[0]*phiUx[1]+((-51192.0*rdxUx[0])-51192.0*rdxLx[0])*phiLy[1]+966.0*rdxLx[0]*phiLx[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*phiC[1]+(58498.28397483125*phiLy[0]-1195.115057222525*phiUx[0])*rdxUx[0]+(1195.115057222525*phiLx[0]-58498.28397483125*phiLy[0])*rdxLx[0])*rdxLyCu[1]+((120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxUxSq[0])*phiUx[1]+((-51552.0*rdxUxSq[0])-287964.0*rdxLx[0]*rdxUx[0]-51552.0*rdxLxSq[0])*phiLy[1]+(120273.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*phiLx[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*phiC[1]+(55172.74642429901*phiUx[0]+60930.08330865796*phiLy[0])*rdxUxSq[0]+(131062.5525579294*phiLx[0]-131062.5525579294*phiUx[0])*rdxLx[0]*rdxUx[0]+((-60930.08330865796*phiLy[0])-55172.74642429901*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+(((-61482.0*rdxUxCu[0])+110061.0*rdxLx[0]*rdxUxSq[0]+122403.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-1944.0*rdxUxCu[0])-132192.0*rdxLx[0]*rdxUxSq[0]-132192.0*rdxLxSq[0]*rdxUx[0]-1944.0*rdxLxCu[0])*phiLy[1]+(122403.0*rdxLx[0]*rdxUxSq[0]+110061.0*rdxLxSq[0]*rdxUx[0]-61482.0*rdxLxCu[0])*phiLx[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*phiC[1]+(57864.35337926103*phiUx[0]+2431.799333826703*phiLy[0])*rdxUxCu[0]+((-136165.1742370272*phiUx[0])+160498.7560325623*phiLy[0]+133234.5442706207*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-133234.5442706207*phiUx[0])-160498.7560325623*phiLy[0]+136165.1742370272*phiLx[0])*rdxLxSq[0]*rdxUx[0]+((-2431.799333826703*phiLy[0])-57864.35337926103*phiLx[0])*rdxLxCu[0])*rdxLy[1]+((-1584.0*rdxUxR4[0])-103032.0*rdxLx[0]*rdxUxCu[0]+205848.0*rdxLxSq[0]*rdxUxSq[0]+3096.0*rdxLxCu[0]*rdxUx[0])*phiUx[1]+(3096.0*rdxLx[0]*rdxUxCu[0]+205848.0*rdxLxSq[0]*rdxUxSq[0]-103032.0*rdxLxCu[0]*rdxUx[0]-1584.0*rdxLxR4[0])*phiLx[1]+(144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0])*phiC[1]+1496.491897739509*phiUx[0]*rdxUxR4[0]+(96897.85037863323*phiUx[0]+3367.106769913895*phiLx[0])*rdxLx[0]*rdxUxCu[0]+(224099.6616864915*phiLx[0]-224099.6616864915*phiUx[0])*rdxLxSq[0]*rdxUxSq[0]+((-3367.106769913895*phiUx[0])-96897.85037863323*phiLx[0])*rdxLxCu[0]*rdxUx[0]-1496.491897739509*phiLx[0]*rdxLxR4[0])*omega-144.0*phiC[1]*rdxUyR4[1]+(((-10086.0*rdxUx[0])-10086.0*rdxLx[0])*phiC[1]-19296.0*phiC[1]*rdxLy[1])*rdxUyCu[1]+((-646704.0*phiC[1]*rdxLySq[1])+((-685848.0*rdxUx[0])-685848.0*rdxLx[0])*phiC[1]*rdxLy[1]+((-19884.0*rdxUxSq[0])-676638.0*rdxLx[0]*rdxUx[0]-19884.0*rdxLxSq[0])*phiC[1])*rdxUySq[1]+((-19296.0*phiC[1]*rdxLyCu[1])+((-685848.0*rdxUx[0])-685848.0*rdxLx[0])*phiC[1]*rdxLySq[1]+((-676638.0*rdxUxSq[0])-1410216.0*rdxLx[0]*rdxUx[0]-676638.0*rdxLxSq[0])*phiC[1]*rdxLy[1]+((-10086.0*rdxUxCu[0])-685848.0*rdxLx[0]*rdxUxSq[0]-685848.0*rdxLxSq[0]*rdxUx[0]-10086.0*rdxLxCu[0])*phiC[1])*rdxUy[1]-144.0*phiC[1]*rdxLyR4[1]+((-10086.0*rdxUx[0])-10086.0*rdxLx[0])*phiC[1]*rdxLyCu[1]+((-19884.0*rdxUxSq[0])-676638.0*rdxLx[0]*rdxUx[0]-19884.0*rdxLxSq[0])*phiC[1]*rdxLySq[1]+((-10086.0*rdxUxCu[0])-685848.0*rdxLx[0]*rdxUxSq[0]-685848.0*rdxLxSq[0]*rdxUx[0]-10086.0*rdxLxCu[0])*phiC[1]*rdxLy[1]+((-144.0*rdxUxR4[0])-19296.0*rdxLx[0]*rdxUxCu[0]-646704.0*rdxLxSq[0]*rdxUxSq[0]-19296.0*rdxLxCu[0]*rdxUx[0]-144.0*rdxLxR4[0])*phiC[1]))/(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0]); 
  phiC[2] = -(1.0*((((24733.68553208356*rdxUx[0]-24733.68553208356*rdxLx[0])*rdxUySq[1]+((4572.614131981835*rdxLx[0]-4572.614131981835*rdxUx[0])*rdxLy[1]+25565.06991971662*rdxUxSq[0]-25565.06991971662*rdxLxSq[0])*rdxUy[1]+(24733.68553208356*rdxUx[0]-24733.68553208356*rdxLx[0])*rdxLySq[1]+(25565.06991971662*rdxUxSq[0]-25565.06991971662*rdxLxSq[0])*rdxLy[1]+831.384387633061*rdxUxCu[0]+54871.36958378201*rdxLx[0]*rdxUxSq[0]-54871.36958378201*rdxLxSq[0]*rdxUx[0]-831.384387633061*rdxLxCu[0])*rho[3]+((-1728.0*rdxUyCu[1])+((-117504.0*rdxLy[1])-69624.0*rdxUx[0]-69624.0*rdxLx[0])*rdxUySq[1]+((-117504.0*rdxLySq[1])+((-251568.0*rdxUx[0])-251568.0*rdxLx[0])*rdxLy[1]-70104.0*rdxUxSq[0]-283728.0*rdxLx[0]*rdxUx[0]-70104.0*rdxLxSq[0])*rdxUy[1]-1728.0*rdxLyCu[1]+((-69624.0*rdxUx[0])-69624.0*rdxLx[0])*rdxLySq[1]+((-70104.0*rdxUxSq[0])-283728.0*rdxLx[0]*rdxUx[0]-70104.0*rdxLxSq[0])*rdxLy[1]-2208.0*rdxUxCu[0]-150144.0*rdxLx[0]*rdxUxSq[0]-150144.0*rdxLxSq[0]*rdxUx[0]-2208.0*rdxLxCu[0])*rho[2]+((63960.0*rdxLx[0]-63960.0*rdxUx[0])*rdxUySq[1]+(63960.0*rdxLxSq[0]-63960.0*rdxUxSq[0])*rdxUy[1]+(63960.0*rdxUx[0]-63960.0*rdxLx[0])*rdxLySq[1]+(63960.0*rdxUxSq[0]-63960.0*rdxLxSq[0])*rdxLy[1])*rho[1]+4323.198815691917*rho[0]*rdxUyCu[1]+(285331.1218356665*rho[0]*rdxLy[1]+(169865.6867998949*rdxUx[0]+169865.6867998949*rdxLx[0])*rho[0])*rdxUySq[1]+((165542.487984203*rdxUxSq[0]+50077.05294843137*rdxLx[0]*rdxUx[0]+165542.487984203*rdxLxSq[0])*rho[0]-285331.1218356665*rho[0]*rdxLySq[1])*rdxUy[1]-4323.198815691917*rho[0]*rdxLyCu[1]+((-169865.6867998949*rdxUx[0])-169865.6867998949*rdxLx[0])*rho[0]*rdxLySq[1]+((-165542.487984203*rdxUxSq[0])-50077.05294843137*rdxLx[0]*rdxUx[0]-165542.487984203*rdxLxSq[0])*rho[0]*rdxLy[1])*omega*volFac+(((47596.75619199274*rdxUx[0]-47596.75619199274*rdxLx[0])*rdxUyCu[1]+((4001.037365484106*rdxUx[0]-4001.037365484106*rdxLx[0])*rdxLy[1]+46869.29485281381*rdxUxSq[0]-46869.29485281381*rdxLxSq[0])*rdxUySq[1]+((90880.70587313897*rdxLx[0]-90880.70587313897*rdxUx[0])*rdxLySq[1]+(91608.1672123179*rdxLxSq[0]-91608.1672123179*rdxUxSq[0])*rdxLy[1]-727.4613391789284*rdxUxCu[0]-48012.44838580926*rdxLx[0]*rdxUxSq[0]+48012.44838580926*rdxLxSq[0]*rdxUx[0]+727.4613391789284*rdxLxCu[0])*rdxUy[1])*phiUy[3]+(1870.614872174387*rdxUx[0]*rdxUyCu[1]+(127201.8113078583*rdxUx[0]*rdxLy[1]+53728.21605078656*rdxUxSq[0]+97012.16573193281*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(127201.8113078583*rdxUx[0]*rdxLySq[1]+(276331.3858395386*rdxUxSq[0]+268329.3111085704*rdxLx[0]*rdxUx[0])*rdxLy[1]+53520.3699538783*rdxUxCu[0]+307144.569706189*rdxLx[0]*rdxUxSq[0]+98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+1870.614872174387*rdxUx[0]*rdxLyCu[1]+(53728.21605078656*rdxUxSq[0]+97012.16573193281*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(53520.3699538783*rdxUxCu[0]+307144.569706189*rdxLx[0]*rdxUxSq[0]+98259.24231338239*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+1662.768775266122*rdxUxR4[0]+114523.1993964541*rdxLx[0]*rdxUxCu[0]+210548.0961680727*rdxLxSq[0]*rdxUxSq[0]+3117.691453623978*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((90880.70587313897*rdxLx[0]-90880.70587313897*rdxUx[0])*rdxLy[1]*rdxUySq[1]+((4001.037365484106*rdxUx[0]-4001.037365484106*rdxLx[0])*rdxLySq[1]+(91608.1672123179*rdxLxSq[0]-91608.1672123179*rdxUxSq[0])*rdxLy[1])*rdxUy[1]+(47596.75619199274*rdxUx[0]-47596.75619199274*rdxLx[0])*rdxLyCu[1]+(46869.29485281381*rdxUxSq[0]-46869.29485281381*rdxLxSq[0])*rdxLySq[1]+((-727.4613391789284*rdxUxCu[0])-48012.44838580926*rdxLx[0]*rdxUxSq[0]+48012.44838580926*rdxLxSq[0]*rdxUx[0]+727.4613391789284*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-1870.614872174387*rdxLx[0]*rdxUyCu[1])+((-127201.8113078583*rdxLx[0]*rdxLy[1])-97012.16573193281*rdxLx[0]*rdxUx[0]-53728.21605078656*rdxLxSq[0])*rdxUySq[1]+((-127201.8113078583*rdxLx[0]*rdxLySq[1])+((-268329.3111085704*rdxLx[0]*rdxUx[0])-276331.3858395386*rdxLxSq[0])*rdxLy[1]-98259.24231338239*rdxLx[0]*rdxUxSq[0]-307144.569706189*rdxLxSq[0]*rdxUx[0]-53520.3699538783*rdxLxCu[0])*rdxUy[1]-1870.614872174387*rdxLx[0]*rdxLyCu[1]+((-97012.16573193281*rdxLx[0]*rdxUx[0])-53728.21605078656*rdxLxSq[0])*rdxLySq[1]+((-98259.24231338239*rdxLx[0]*rdxUxSq[0])-307144.569706189*rdxLxSq[0]*rdxUx[0]-53520.3699538783*rdxLxCu[0])*rdxLy[1]-3117.691453623978*rdxLx[0]*rdxUxCu[0]-210548.0961680727*rdxLxSq[0]*rdxUxSq[0]-114523.1993964541*rdxLxCu[0]*rdxUx[0]-1662.768775266122*rdxLxR4[0])*phiLx[3]+((-1584.0*rdxUyR4[1])+((-103032.0*rdxLy[1])-61482.0*rdxUx[0]-61482.0*rdxLx[0])*rdxUyCu[1]+(205848.0*rdxLySq[1]+(110061.0*rdxUx[0]+110061.0*rdxLx[0])*rdxLy[1]-58932.0*rdxUxSq[0]+97026.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*rdxUySq[1]+(3096.0*rdxLyCu[1]+(122403.0*rdxUx[0]+122403.0*rdxLx[0])*rdxLySq[1]+(120273.0*rdxUxSq[0]+151236.0*rdxLx[0]*rdxUx[0]+120273.0*rdxLxSq[0])*rdxLy[1]+966.0*rdxUxCu[0]+65688.0*rdxLx[0]*rdxUxSq[0]+65688.0*rdxLxSq[0]*rdxUx[0]+966.0*rdxLxCu[0])*rdxUy[1])*phiUy[2]+((-1944.0*rdxUx[0]*rdxUyCu[1])+((-132192.0*rdxUx[0]*rdxLy[1])-51552.0*rdxUxSq[0]-105102.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-132192.0*rdxUx[0]*rdxLySq[1])+((-287964.0*rdxUxSq[0])-278064.0*rdxLx[0]*rdxUx[0])*rdxLy[1]-51192.0*rdxUxCu[0]-319194.0*rdxLx[0]*rdxUxSq[0]-106542.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-1944.0*rdxUx[0]*rdxLyCu[1]+((-51552.0*rdxUxSq[0])-105102.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-51192.0*rdxUxCu[0])-319194.0*rdxLx[0]*rdxUxSq[0]-106542.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-1584.0*rdxUxR4[0]-109512.0*rdxLx[0]*rdxUxCu[0]-228312.0*rdxLxSq[0]*rdxUxSq[0]-3384.0*rdxLxCu[0]*rdxUx[0])*phiUx[2]+(3096.0*rdxLy[1]*rdxUyCu[1]+(205848.0*rdxLySq[1]+(122403.0*rdxUx[0]+122403.0*rdxLx[0])*rdxLy[1])*rdxUySq[1]+((-103032.0*rdxLyCu[1])+(110061.0*rdxUx[0]+110061.0*rdxLx[0])*rdxLySq[1]+(120273.0*rdxUxSq[0]+151236.0*rdxLx[0]*rdxUx[0]+120273.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]-1584.0*rdxLyR4[1]+((-61482.0*rdxUx[0])-61482.0*rdxLx[0])*rdxLyCu[1]+((-58932.0*rdxUxSq[0])+97026.0*rdxLx[0]*rdxUx[0]-58932.0*rdxLxSq[0])*rdxLySq[1]+(966.0*rdxUxCu[0]+65688.0*rdxLx[0]*rdxUxSq[0]+65688.0*rdxLxSq[0]*rdxUx[0]+966.0*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-1944.0*rdxLx[0]*rdxUyCu[1])+((-132192.0*rdxLx[0]*rdxLy[1])-105102.0*rdxLx[0]*rdxUx[0]-51552.0*rdxLxSq[0])*rdxUySq[1]+((-132192.0*rdxLx[0]*rdxLySq[1])+((-278064.0*rdxLx[0]*rdxUx[0])-287964.0*rdxLxSq[0])*rdxLy[1]-106542.0*rdxLx[0]*rdxUxSq[0]-319194.0*rdxLxSq[0]*rdxUx[0]-51192.0*rdxLxCu[0])*rdxUy[1]-1944.0*rdxLx[0]*rdxLyCu[1]+((-105102.0*rdxLx[0]*rdxUx[0])-51552.0*rdxLxSq[0])*rdxLySq[1]+((-106542.0*rdxLx[0]*rdxUxSq[0])-319194.0*rdxLxSq[0]*rdxUx[0]-51192.0*rdxLxCu[0])*rdxLy[1]-3384.0*rdxLx[0]*rdxUxCu[0]-228312.0*rdxLxSq[0]*rdxUxSq[0]-109512.0*rdxLxCu[0]*rdxUx[0]-1584.0*rdxLxR4[0])*phiLx[2]+(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0])*phiC[2]+1496.491897739509*phiUy[0]*rdxUyR4[1]+((96897.85037863323*phiUy[0]+3367.106769913895*phiLy[0])*rdxLy[1]+(45180.0*rdxLx[0]-45180.0*rdxUx[0])*phiUy[1]-2340.0*rdxUx[0]*phiUx[1]+2340.0*rdxLx[0]*phiLx[1]+(57864.35337926103*phiUy[0]+2431.799333826703*phiUx[0])*rdxUx[0]+(57864.35337926103*phiUy[0]+2431.799333826703*phiLx[0])*rdxLx[0])*rdxUyCu[1]+((224099.6616864915*phiLy[0]-224099.6616864915*phiUy[0])*rdxLySq[1]+((4950.0*rdxLx[0]-4950.0*rdxUx[0])*phiUy[1]-154440.0*rdxUx[0]*phiUx[1]+(98730.0*rdxLx[0]-98730.0*rdxUx[0])*phiLy[1]+154440.0*rdxLx[0]*phiLx[1]+((-136165.1742370272*phiUy[0])+160498.7560325623*phiUx[0]+133234.5442706207*phiLy[0])*rdxUx[0]+((-136165.1742370272*phiUy[0])+133234.5442706207*phiLy[0]+160498.7560325623*phiLx[0])*rdxLx[0])*rdxLy[1]+(44280.0*rdxLxSq[0]-44280.0*rdxUxSq[0])*phiUy[1]+((-63960.0*rdxUxSq[0])-119925.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(119925.0*rdxLx[0]*rdxUx[0]+63960.0*rdxLxSq[0])*phiLx[1]+(55172.74642429901*phiUy[0]+60930.08330865796*phiUx[0])*rdxUxSq[0]+((-125403.9425696018*phiUy[0])+130168.8143412238*phiUx[0]+130168.8143412238*phiLx[0])*rdxLx[0]*rdxUx[0]+(55172.74642429901*phiUy[0]+60930.08330865796*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+(((-3367.106769913895*phiUy[0])-96897.85037863323*phiLy[0])*rdxLyCu[1]+((98730.0*rdxUx[0]-98730.0*rdxLx[0])*phiUy[1]+154440.0*rdxUx[0]*phiUx[1]+(4950.0*rdxUx[0]-4950.0*rdxLx[0])*phiLy[1]-154440.0*rdxLx[0]*phiLx[1]+((-133234.5442706207*phiUy[0])-160498.7560325623*phiUx[0]+136165.1742370272*phiLy[0])*rdxUx[0]+((-133234.5442706207*phiUy[0])+136165.1742370272*phiLy[0]-160498.7560325623*phiLx[0])*rdxLx[0])*rdxLySq[1]+((99630.0*rdxUxSq[0]-99630.0*rdxLxSq[0])*phiUy[1]+(99630.0*rdxLxSq[0]-99630.0*rdxUxSq[0])*phiLy[1]+(131062.5525579294*phiLy[0]-131062.5525579294*phiUy[0])*rdxUxSq[0]+(181740.6271365871*phiLy[0]-181740.6271365871*phiUy[0])*rdxLx[0]*rdxUx[0]+(131062.5525579294*phiLy[0]-131062.5525579294*phiUy[0])*rdxLxSq[0])*rdxLy[1]+(900.0*rdxUxCu[0]+59400.0*rdxLx[0]*rdxUxSq[0]-59400.0*rdxLxSq[0]*rdxUx[0]-900.0*rdxLxCu[0])*phiUy[1]+((-61620.0*rdxUxCu[0])-27105.0*rdxLx[0]*rdxUxSq[0]-117585.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(117585.0*rdxLx[0]*rdxUxSq[0]+27105.0*rdxLxSq[0]*rdxUx[0]+61620.0*rdxLxCu[0])*phiLx[1]+(58498.28397483125*phiUx[0]-1195.115057222525*phiUy[0])*rdxUxCu[0]+((-81267.82389113172*phiUy[0])+28168.34228349264*phiUx[0]+127737.0150073971*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-81267.82389113172*phiUy[0])+127737.0150073971*phiUx[0]+28168.34228349264*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(58498.28397483125*phiLx[0]-1195.115057222525*phiUy[0])*rdxLxCu[0])*rdxUy[1]-1496.491897739509*phiLy[0]*rdxLyR4[1]+(2340.0*rdxUx[0]*phiUx[1]+(45180.0*rdxUx[0]-45180.0*rdxLx[0])*phiLy[1]-2340.0*rdxLx[0]*phiLx[1]+((-2431.799333826703*phiUx[0])-57864.35337926103*phiLy[0])*rdxUx[0]+((-57864.35337926103*phiLy[0])-2431.799333826703*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((63960.0*rdxUxSq[0]+119925.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(44280.0*rdxUxSq[0]-44280.0*rdxLxSq[0])*phiLy[1]+((-119925.0*rdxLx[0]*rdxUx[0])-63960.0*rdxLxSq[0])*phiLx[1]+((-60930.08330865796*phiUx[0])-55172.74642429901*phiLy[0])*rdxUxSq[0]+((-130168.8143412238*phiUx[0])+125403.9425696018*phiLy[0]-130168.8143412238*phiLx[0])*rdxLx[0]*rdxUx[0]+((-55172.74642429901*phiLy[0])-60930.08330865796*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((61620.0*rdxUxCu[0]+27105.0*rdxLx[0]*rdxUxSq[0]+117585.0*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-900.0*rdxUxCu[0])-59400.0*rdxLx[0]*rdxUxSq[0]+59400.0*rdxLxSq[0]*rdxUx[0]+900.0*rdxLxCu[0])*phiLy[1]+((-117585.0*rdxLx[0]*rdxUxSq[0])-27105.0*rdxLxSq[0]*rdxUx[0]-61620.0*rdxLxCu[0])*phiLx[1]+(1195.115057222525*phiLy[0]-58498.28397483125*phiUx[0])*rdxUxCu[0]+((-28168.34228349264*phiUx[0])+81267.82389113172*phiLy[0]-127737.0150073971*phiLx[0])*rdxLx[0]*rdxUxSq[0]+((-127737.0150073971*phiUx[0])+81267.82389113172*phiLy[0]-28168.34228349264*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(1195.115057222525*phiLy[0]-58498.28397483125*phiLx[0])*rdxLxCu[0])*rdxLy[1])*omega+((-144.0*rdxUyR4[1])+((-19296.0*rdxLy[1])-10086.0*rdxUx[0]-10086.0*rdxLx[0])*rdxUyCu[1]+((-646704.0*rdxLySq[1])+((-685848.0*rdxUx[0])-685848.0*rdxLx[0])*rdxLy[1]-19884.0*rdxUxSq[0]-676638.0*rdxLx[0]*rdxUx[0]-19884.0*rdxLxSq[0])*rdxUySq[1]+((-19296.0*rdxLyCu[1])+((-685848.0*rdxUx[0])-685848.0*rdxLx[0])*rdxLySq[1]+((-676638.0*rdxUxSq[0])-1410216.0*rdxLx[0]*rdxUx[0]-676638.0*rdxLxSq[0])*rdxLy[1]-10086.0*rdxUxCu[0]-685848.0*rdxLx[0]*rdxUxSq[0]-685848.0*rdxLxSq[0]*rdxUx[0]-10086.0*rdxLxCu[0])*rdxUy[1]-144.0*rdxLyR4[1]+((-10086.0*rdxUx[0])-10086.0*rdxLx[0])*rdxLyCu[1]+((-19884.0*rdxUxSq[0])-676638.0*rdxLx[0]*rdxUx[0]-19884.0*rdxLxSq[0])*rdxLySq[1]+((-10086.0*rdxUxCu[0])-685848.0*rdxLx[0]*rdxUxSq[0]-685848.0*rdxLxSq[0]*rdxUx[0]-10086.0*rdxLxCu[0])*rdxLy[1]-144.0*rdxUxR4[0]-19296.0*rdxLx[0]*rdxUxCu[0]-646704.0*rdxLxSq[0]*rdxUxSq[0]-19296.0*rdxLxCu[0]*rdxUx[0]-144.0*rdxLxR4[0])*phiC[2]))/(144.0*rdxUyR4[1]+(19296.0*rdxLy[1]+10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxUyCu[1]+(646704.0*rdxLySq[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLy[1]+19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxUySq[1]+(19296.0*rdxLyCu[1]+(685848.0*rdxUx[0]+685848.0*rdxLx[0])*rdxLySq[1]+(676638.0*rdxUxSq[0]+1410216.0*rdxLx[0]*rdxUx[0]+676638.0*rdxLxSq[0])*rdxLy[1]+10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxUy[1]+144.0*rdxLyR4[1]+(10086.0*rdxUx[0]+10086.0*rdxLx[0])*rdxLyCu[1]+(19884.0*rdxUxSq[0]+676638.0*rdxLx[0]*rdxUx[0]+19884.0*rdxLxSq[0])*rdxLySq[1]+(10086.0*rdxUxCu[0]+685848.0*rdxLx[0]*rdxUxSq[0]+685848.0*rdxLxSq[0]*rdxUx[0]+10086.0*rdxLxCu[0])*rdxLy[1]+144.0*rdxUxR4[0]+19296.0*rdxLx[0]*rdxUxCu[0]+646704.0*rdxLxSq[0]*rdxUxSq[0]+19296.0*rdxLxCu[0]*rdxUx[0]+144.0*rdxLxR4[0]); 
  phiC[3] = (((288.0*rdxUyCu[1]+(19584.0*rdxLy[1]+9224.0*rdxUx[0]+9224.0*rdxLx[0])*rdxUySq[1]+(19584.0*rdxLySq[1]+(42368.0*rdxUx[0]+42368.0*rdxLx[0])*rdxLy[1]+9224.0*rdxUxSq[0]+42368.0*rdxLx[0]*rdxUx[0]+9224.0*rdxLxSq[0])*rdxUy[1]+288.0*rdxLyCu[1]+(9224.0*rdxUx[0]+9224.0*rdxLx[0])*rdxLySq[1]+(9224.0*rdxUxSq[0]+42368.0*rdxLx[0]*rdxUx[0]+9224.0*rdxLxSq[0])*rdxLy[1]+288.0*rdxUxCu[0]+19584.0*rdxLx[0]*rdxUxSq[0]+19584.0*rdxLxSq[0]*rdxUx[0]+288.0*rdxLxCu[0])*rho[3]+((21435.86079447242*rdxLx[0]-21435.86079447242*rdxUx[0])*rdxUySq[1]+((3962.932247717591*rdxUx[0]-3962.932247717591*rdxLx[0])*rdxLy[1]-22156.39393042108*rdxUxSq[0]+22156.39393042108*rdxLxSq[0])*rdxUy[1]+(21435.86079447242*rdxLx[0]-21435.86079447242*rdxUx[0])*rdxLySq[1]+(22156.39393042108*rdxLxSq[0]-22156.39393042108*rdxUxSq[0])*rdxLy[1]-720.5331359486529*rdxUxCu[0]-47555.18697261109*rdxLx[0]*rdxUxSq[0]+47555.18697261109*rdxLxSq[0]*rdxUx[0]+720.5331359486529*rdxLxCu[0])*rho[2]+((-720.5331359486529*rdxUyCu[1])+((-47555.18697261109*rdxLy[1])-22156.39393042108*rdxUx[0]-22156.39393042108*rdxLx[0])*rdxUySq[1]+(47555.18697261109*rdxLySq[1]-21435.86079447242*rdxUxSq[0]+3962.932247717591*rdxLx[0]*rdxUx[0]-21435.86079447242*rdxLxSq[0])*rdxUy[1]+720.5331359486529*rdxLyCu[1]+(22156.39393042108*rdxUx[0]+22156.39393042108*rdxLx[0])*rdxLySq[1]+(21435.86079447242*rdxUxSq[0]-3962.932247717591*rdxLx[0]*rdxUx[0]+21435.86079447242*rdxLxSq[0])*rdxLy[1])*rho[1]+(55432.0*rdxUx[0]-55432.0*rdxLx[0])*rho[0]*rdxUySq[1]+(55432.0*rdxUxSq[0]-55432.0*rdxLxSq[0])*rho[0]*rdxUy[1]+(55432.0*rdxLx[0]-55432.0*rdxUx[0])*rho[0]*rdxLySq[1]+(55432.0*rdxLxSq[0]-55432.0*rdxUxSq[0])*rho[0]*rdxLy[1])*omega*volFac+((528.0*rdxUyR4[1]+(34344.0*rdxLy[1]+15914.0*rdxUx[0]+15914.0*rdxLx[0])*rdxUyCu[1]+((-68616.0*rdxLySq[1])+((-37072.0*rdxUx[0])-37072.0*rdxLx[0])*rdxLy[1]+15134.0*rdxUxSq[0]-41362.0*rdxLx[0]*rdxUx[0]+15134.0*rdxLxSq[0])*rdxUySq[1]+((-1032.0*rdxLyCu[1])+((-32056.0*rdxUx[0])-32056.0*rdxLx[0])*rdxLySq[1]+((-31276.0*rdxUxSq[0])-32782.0*rdxLx[0]*rdxUx[0]-31276.0*rdxLxSq[0])*rdxLy[1]-252.0*rdxUxCu[0]-17136.0*rdxLx[0]*rdxUxSq[0]-17136.0*rdxLxSq[0]*rdxUx[0]-252.0*rdxLxCu[0])*rdxUy[1])*phiUy[3]+((-252.0*rdxUx[0]*rdxUyCu[1])+((-17136.0*rdxUx[0]*rdxLy[1])+15134.0*rdxUxSq[0]-31276.0*rdxLx[0]*rdxUx[0])*rdxUySq[1]+((-17136.0*rdxUx[0]*rdxLySq[1])+((-41362.0*rdxUxSq[0])-32782.0*rdxLx[0]*rdxUx[0])*rdxLy[1]+15914.0*rdxUxCu[0]-37072.0*rdxLx[0]*rdxUxSq[0]-32056.0*rdxLxSq[0]*rdxUx[0])*rdxUy[1]-252.0*rdxUx[0]*rdxLyCu[1]+(15134.0*rdxUxSq[0]-31276.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+(15914.0*rdxUxCu[0]-37072.0*rdxLx[0]*rdxUxSq[0]-32056.0*rdxLxSq[0]*rdxUx[0])*rdxLy[1]+528.0*rdxUxR4[0]+34344.0*rdxLx[0]*rdxUxCu[0]-68616.0*rdxLxSq[0]*rdxUxSq[0]-1032.0*rdxLxCu[0]*rdxUx[0])*phiUx[3]+((-1032.0*rdxLy[1]*rdxUyCu[1])+(((-32056.0*rdxUx[0])-32056.0*rdxLx[0])*rdxLy[1]-68616.0*rdxLySq[1])*rdxUySq[1]+(34344.0*rdxLyCu[1]+((-37072.0*rdxUx[0])-37072.0*rdxLx[0])*rdxLySq[1]+((-31276.0*rdxUxSq[0])-32782.0*rdxLx[0]*rdxUx[0]-31276.0*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+528.0*rdxLyR4[1]+(15914.0*rdxUx[0]+15914.0*rdxLx[0])*rdxLyCu[1]+(15134.0*rdxUxSq[0]-41362.0*rdxLx[0]*rdxUx[0]+15134.0*rdxLxSq[0])*rdxLySq[1]+((-252.0*rdxUxCu[0])-17136.0*rdxLx[0]*rdxUxSq[0]-17136.0*rdxLxSq[0]*rdxUx[0]-252.0*rdxLxCu[0])*rdxLy[1])*phiLy[3]+((-252.0*rdxLx[0]*rdxUyCu[1])+((-17136.0*rdxLx[0]*rdxLy[1])-31276.0*rdxLx[0]*rdxUx[0]+15134.0*rdxLxSq[0])*rdxUySq[1]+((-17136.0*rdxLx[0]*rdxLySq[1])+((-32782.0*rdxLx[0]*rdxUx[0])-41362.0*rdxLxSq[0])*rdxLy[1]-32056.0*rdxLx[0]*rdxUxSq[0]-37072.0*rdxLxSq[0]*rdxUx[0]+15914.0*rdxLxCu[0])*rdxUy[1]-252.0*rdxLx[0]*rdxLyCu[1]+(15134.0*rdxLxSq[0]-31276.0*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-32056.0*rdxLx[0]*rdxUxSq[0])-37072.0*rdxLxSq[0]*rdxUx[0]+15914.0*rdxLxCu[0])*rdxLy[1]-1032.0*rdxLx[0]*rdxUxCu[0]-68616.0*rdxLxSq[0]*rdxUxSq[0]+34344.0*rdxLxCu[0]*rdxUx[0]+528.0*rdxLxR4[0])*phiLx[3]+((-48.0*rdxUyR4[1])+((-6432.0*rdxLy[1])-3362.0*rdxUx[0]-3362.0*rdxLx[0])*rdxUyCu[1]+((-215568.0*rdxLySq[1])+((-228616.0*rdxUx[0])-228616.0*rdxLx[0])*rdxLy[1]-6628.0*rdxUxSq[0]-225546.0*rdxLx[0]*rdxUx[0]-6628.0*rdxLxSq[0])*rdxUySq[1]+((-6432.0*rdxLyCu[1])+((-228616.0*rdxUx[0])-228616.0*rdxLx[0])*rdxLySq[1]+((-225546.0*rdxUxSq[0])-470072.0*rdxLx[0]*rdxUx[0]-225546.0*rdxLxSq[0])*rdxLy[1]-3362.0*rdxUxCu[0]-228616.0*rdxLx[0]*rdxUxSq[0]-228616.0*rdxLxSq[0]*rdxUx[0]-3362.0*rdxLxCu[0])*rdxUy[1]-48.0*rdxLyR4[1]+((-3362.0*rdxUx[0])-3362.0*rdxLx[0])*rdxLyCu[1]+((-6628.0*rdxUxSq[0])-225546.0*rdxLx[0]*rdxUx[0]-6628.0*rdxLxSq[0])*rdxLySq[1]+((-3362.0*rdxUxCu[0])-228616.0*rdxLx[0]*rdxUxSq[0]-228616.0*rdxLxSq[0]*rdxUx[0]-3362.0*rdxLxCu[0])*rdxLy[1]-48.0*rdxUxR4[0]-6432.0*rdxLx[0]*rdxUxCu[0]-215568.0*rdxLxSq[0]*rdxUxSq[0]-6432.0*rdxLxCu[0]*rdxUx[0]-48.0*rdxLxR4[0])*phiC[3]+((20625.26101653019*rdxLx[0]-20625.26101653019*rdxUx[0])*rdxUyCu[1]+((1733.782858376446*rdxLx[0]-1733.782858376446*rdxUx[0])*rdxLy[1]-20310.02776955265*rdxUxSq[0]+20310.02776955265*rdxLxSq[0])*rdxUySq[1]+((39381.63921169356*rdxUx[0]-39381.63921169356*rdxLx[0])*rdxLySq[1]+(39696.8724586711*rdxUxSq[0]-39696.8724586711*rdxLxSq[0])*rdxLy[1]+315.2332469775357*rdxUxCu[0]+20805.39430051735*rdxLx[0]*rdxUxSq[0]-20805.39430051735*rdxLxSq[0]*rdxUx[0]-315.2332469775357*rdxLxCu[0])*rdxUy[1])*phiUy[2]+(311.7691453623978*rdxUx[0]*rdxUyCu[1]+(21200.30188464305*rdxUx[0]*rdxLy[1]-14130.0704881469*rdxUxSq[0]+34100.61629941605*rdxLx[0]*rdxUx[0])*rdxUySq[1]+(21200.30188464305*rdxUx[0]*rdxLySq[1]+(50323.00416310616*rdxUxSq[0]+41406.40660574158*rdxLx[0]*rdxUx[0])*rdxLy[1]-14940.67026608914*rdxUxCu[0]+45864.70538442387*rdxLx[0]*rdxUxSq[0]+34911.21607735829*rdxLxSq[0]*rdxUx[0])*rdxUy[1]+311.7691453623978*rdxUx[0]*rdxLyCu[1]+(34100.61629941605*rdxLx[0]*rdxUx[0]-14130.0704881469*rdxUxSq[0])*rdxLySq[1]+((-14940.67026608914*rdxUxCu[0])+45864.70538442387*rdxLx[0]*rdxUxSq[0]+34911.21607735829*rdxLxSq[0]*rdxUx[0])*rdxLy[1]-498.8306325798365*rdxUxR4[0]-32299.28345954441*rdxLx[0]*rdxUxCu[0]+74699.88722883051*rdxLxSq[0]*rdxUxSq[0]+1122.368923304632*rdxLxCu[0]*rdxUx[0])*phiUx[2]+((39381.63921169356*rdxUx[0]-39381.63921169356*rdxLx[0])*rdxLy[1]*rdxUySq[1]+((1733.782858376446*rdxLx[0]-1733.782858376446*rdxUx[0])*rdxLySq[1]+(39696.8724586711*rdxUxSq[0]-39696.8724586711*rdxLxSq[0])*rdxLy[1])*rdxUy[1]+(20625.26101653019*rdxLx[0]-20625.26101653019*rdxUx[0])*rdxLyCu[1]+(20310.02776955265*rdxLxSq[0]-20310.02776955265*rdxUxSq[0])*rdxLySq[1]+(315.2332469775357*rdxUxCu[0]+20805.39430051735*rdxLx[0]*rdxUxSq[0]-20805.39430051735*rdxLxSq[0]*rdxUx[0]-315.2332469775357*rdxLxCu[0])*rdxLy[1])*phiLy[2]+((-311.7691453623978*rdxLx[0]*rdxUyCu[1])+((-21200.30188464305*rdxLx[0]*rdxLy[1])-34100.61629941605*rdxLx[0]*rdxUx[0]+14130.0704881469*rdxLxSq[0])*rdxUySq[1]+((-21200.30188464305*rdxLx[0]*rdxLySq[1])+((-41406.40660574158*rdxLx[0]*rdxUx[0])-50323.00416310616*rdxLxSq[0])*rdxLy[1]-34911.21607735829*rdxLx[0]*rdxUxSq[0]-45864.70538442387*rdxLxSq[0]*rdxUx[0]+14940.67026608914*rdxLxCu[0])*rdxUy[1]-311.7691453623978*rdxLx[0]*rdxLyCu[1]+(14130.0704881469*rdxLxSq[0]-34100.61629941605*rdxLx[0]*rdxUx[0])*rdxLySq[1]+((-34911.21607735829*rdxLx[0]*rdxUxSq[0])-45864.70538442387*rdxLxSq[0]*rdxUx[0]+14940.67026608914*rdxLxCu[0])*rdxLy[1]-1122.368923304632*rdxLx[0]*rdxUxCu[0]-74699.88722883051*rdxLxSq[0]*rdxUxSq[0]+32299.28345954441*rdxLxCu[0]*rdxUx[0]+498.8306325798365*rdxLxR4[0])*phiLx[2]-498.8306325798365*phiUy[1]*rdxUyR4[1]+(((-32299.28345954441*phiUy[1])-1122.368923304632*phiLy[1])*rdxLy[1]+((-14940.67026608914*rdxUx[0])-14940.67026608914*rdxLx[0])*phiUy[1]+315.2332469775357*rdxUx[0]*phiUx[1]+315.2332469775357*rdxLx[0]*phiLx[1]+(19578.0*phiUy[0]-390.0*phiUx[0])*rdxUx[0]+(390.0*phiLx[0]-19578.0*phiUy[0])*rdxLx[0])*rdxUyCu[1]+((74699.88722883051*phiUy[1]-74699.88722883051*phiLy[1])*rdxLySq[1]+((45864.70538442387*rdxUx[0]+45864.70538442387*rdxLx[0])*phiUy[1]+20805.39430051735*rdxUx[0]*phiUx[1]+((-34911.21607735829*rdxUx[0])-34911.21607735829*rdxLx[0])*phiLy[1]+20805.39430051735*rdxLx[0]*phiLx[1]+(2145.0*phiUy[0]-25740.0*phiUx[0]+42783.0*phiLy[0])*rdxUx[0]+((-2145.0*phiUy[0])-42783.0*phiLy[0]+25740.0*phiLx[0])*rdxLx[0])*rdxLy[1]+((-14130.0704881469*rdxUxSq[0])+50323.00416310616*rdxLx[0]*rdxUx[0]-14130.0704881469*rdxLxSq[0])*phiUy[1]+(39696.8724586711*rdxLx[0]*rdxUx[0]-20310.02776955265*rdxUxSq[0])*phiUx[1]+(39696.8724586711*rdxLx[0]*rdxUx[0]-20310.02776955265*rdxLxSq[0])*phiLx[1]+(19188.0*phiUy[0]+19188.0*phiUx[0])*rdxUxSq[0]+(43173.0*phiLx[0]-43173.0*phiUx[0])*rdxLx[0]*rdxUx[0]+((-19188.0*phiUy[0])-19188.0*phiLx[0])*rdxLxSq[0])*rdxUySq[1]+((1122.368923304632*phiUy[1]+32299.28345954441*phiLy[1])*rdxLyCu[1]+((34911.21607735829*rdxUx[0]+34911.21607735829*rdxLx[0])*phiUy[1]-20805.39430051735*rdxUx[0]*phiUx[1]+((-45864.70538442387*rdxUx[0])-45864.70538442387*rdxLx[0])*phiLy[1]-20805.39430051735*rdxLx[0]*phiLx[1]+((-42783.0*phiUy[0])+25740.0*phiUx[0]-2145.0*phiLy[0])*rdxUx[0]+(42783.0*phiUy[0]+2145.0*phiLy[0]-25740.0*phiLx[0])*rdxLx[0])*rdxLySq[1]+((34100.61629941605*rdxUxSq[0]+41406.40660574158*rdxLx[0]*rdxUx[0]+34100.61629941605*rdxLxSq[0])*phiUy[1]+((-34100.61629941605*rdxUxSq[0])-41406.40660574158*rdxLx[0]*rdxUx[0]-34100.61629941605*rdxLxSq[0])*phiLy[1]+(43173.0*phiLy[0]-43173.0*phiUy[0])*rdxUxSq[0]+(43173.0*phiUy[0]-43173.0*phiLy[0])*rdxLxSq[0])*rdxLy[1]+(311.7691453623978*rdxUxCu[0]+21200.30188464305*rdxLx[0]*rdxUxSq[0]+21200.30188464305*rdxLxSq[0]*rdxUx[0]+311.7691453623978*rdxLxCu[0])*phiUy[1]+((-20625.26101653019*rdxUxCu[0])-1733.782858376446*rdxLx[0]*rdxUxSq[0]+39381.63921169356*rdxLxSq[0]*rdxUx[0])*phiUx[1]+(39381.63921169356*rdxLx[0]*rdxUxSq[0]-1733.782858376446*rdxLxSq[0]*rdxUx[0]-20625.26101653019*rdxLxCu[0])*phiLx[1]+(19578.0*phiUx[0]-390.0*phiUy[0])*rdxUxCu[0]+((-25740.0*phiUy[0])+2145.0*phiUx[0]+42783.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(25740.0*phiUy[0]-42783.0*phiUx[0]-2145.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(390.0*phiUy[0]-19578.0*phiLx[0])*rdxLxCu[0])*rdxUy[1]+498.8306325798365*phiLy[1]*rdxLyR4[1]+((-315.2332469775357*rdxUx[0]*phiUx[1])+(14940.67026608914*rdxUx[0]+14940.67026608914*rdxLx[0])*phiLy[1]-315.2332469775357*rdxLx[0]*phiLx[1]+(390.0*phiUx[0]-19578.0*phiLy[0])*rdxUx[0]+(19578.0*phiLy[0]-390.0*phiLx[0])*rdxLx[0])*rdxLyCu[1]+((20310.02776955265*rdxUxSq[0]-39696.8724586711*rdxLx[0]*rdxUx[0])*phiUx[1]+(14130.0704881469*rdxUxSq[0]-50323.00416310616*rdxLx[0]*rdxUx[0]+14130.0704881469*rdxLxSq[0])*phiLy[1]+(20310.02776955265*rdxLxSq[0]-39696.8724586711*rdxLx[0]*rdxUx[0])*phiLx[1]+((-19188.0*phiUx[0])-19188.0*phiLy[0])*rdxUxSq[0]+(43173.0*phiUx[0]-43173.0*phiLx[0])*rdxLx[0]*rdxUx[0]+(19188.0*phiLy[0]+19188.0*phiLx[0])*rdxLxSq[0])*rdxLySq[1]+((20625.26101653019*rdxUxCu[0]+1733.782858376446*rdxLx[0]*rdxUxSq[0]-39381.63921169356*rdxLxSq[0]*rdxUx[0])*phiUx[1]+((-311.7691453623978*rdxUxCu[0])-21200.30188464305*rdxLx[0]*rdxUxSq[0]-21200.30188464305*rdxLxSq[0]*rdxUx[0]-311.7691453623978*rdxLxCu[0])*phiLy[1]+((-39381.63921169356*rdxLx[0]*rdxUxSq[0])+1733.782858376446*rdxLxSq[0]*rdxUx[0]+20625.26101653019*rdxLxCu[0])*phiLx[1]+(390.0*phiLy[0]-19578.0*phiUx[0])*rdxUxCu[0]+((-2145.0*phiUx[0])+25740.0*phiLy[0]-42783.0*phiLx[0])*rdxLx[0]*rdxUxSq[0]+(42783.0*phiUx[0]-25740.0*phiLy[0]+2145.0*phiLx[0])*rdxLxSq[0]*rdxUx[0]+(19578.0*phiLx[0]-390.0*phiLy[0])*rdxLxCu[0])*rdxLy[1])*omega+(48.0*rdxUyR4[1]+(6432.0*rdxLy[1]+3362.0*rdxUx[0]+3362.0*rdxLx[0])*rdxUyCu[1]+(215568.0*rdxLySq[1]+(228616.0*rdxUx[0]+228616.0*rdxLx[0])*rdxLy[1]+6628.0*rdxUxSq[0]+225546.0*rdxLx[0]*rdxUx[0]+6628.0*rdxLxSq[0])*rdxUySq[1]+(6432.0*rdxLyCu[1]+(228616.0*rdxUx[0]+228616.0*rdxLx[0])*rdxLySq[1]+(225546.0*rdxUxSq[0]+470072.0*rdxLx[0]*rdxUx[0]+225546.0*rdxLxSq[0])*rdxLy[1]+3362.0*rdxUxCu[0]+228616.0*rdxLx[0]*rdxUxSq[0]+228616.0*rdxLxSq[0]*rdxUx[0]+3362.0*rdxLxCu[0])*rdxUy[1]+48.0*rdxLyR4[1]+(3362.0*rdxUx[0]+3362.0*rdxLx[0])*rdxLyCu[1]+(6628.0*rdxUxSq[0]+225546.0*rdxLx[0]*rdxUx[0]+6628.0*rdxLxSq[0])*rdxLySq[1]+(3362.0*rdxUxCu[0]+228616.0*rdxLx[0]*rdxUxSq[0]+228616.0*rdxLxSq[0]*rdxUx[0]+3362.0*rdxLxCu[0])*rdxLy[1]+48.0*rdxUxR4[0]+6432.0*rdxLx[0]*rdxUxCu[0]+215568.0*rdxLxSq[0]*rdxUxSq[0]+6432.0*rdxLxCu[0]*rdxUx[0]+48.0*rdxLxR4[0])*phiC[3])/(48.0*rdxUyR4[1]+(6432.0*rdxLy[1]+3362.0*rdxUx[0]+3362.0*rdxLx[0])*rdxUyCu[1]+(215568.0*rdxLySq[1]+(228616.0*rdxUx[0]+228616.0*rdxLx[0])*rdxLy[1]+6628.0*rdxUxSq[0]+225546.0*rdxLx[0]*rdxUx[0]+6628.0*rdxLxSq[0])*rdxUySq[1]+(6432.0*rdxLyCu[1]+(228616.0*rdxUx[0]+228616.0*rdxLx[0])*rdxLySq[1]+(225546.0*rdxUxSq[0]+470072.0*rdxLx[0]*rdxUx[0]+225546.0*rdxLxSq[0])*rdxLy[1]+3362.0*rdxUxCu[0]+228616.0*rdxLx[0]*rdxUxSq[0]+228616.0*rdxLxSq[0]*rdxUx[0]+3362.0*rdxLxCu[0])*rdxUy[1]+48.0*rdxLyR4[1]+(3362.0*rdxUx[0]+3362.0*rdxLx[0])*rdxLyCu[1]+(6628.0*rdxUxSq[0]+225546.0*rdxLx[0]*rdxUx[0]+6628.0*rdxLxSq[0])*rdxLySq[1]+(3362.0*rdxUxCu[0]+228616.0*rdxLx[0]*rdxUxSq[0]+228616.0*rdxLxSq[0]*rdxUx[0]+3362.0*rdxLxCu[0])*rdxLy[1]+48.0*rdxUxR4[0]+6432.0*rdxLx[0]*rdxUxCu[0]+215568.0*rdxLxSq[0]*rdxUxSq[0]+6432.0*rdxLxCu[0]*rdxUx[0]+48.0*rdxLxR4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((1.732050807568877*rdxCp2[0]*rho[1]+6.0*rho[0]*rdxCp2[1]+50.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+5.196152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]+9.0*phiUx[0]+75.0*phiLy[0]-177.0*phiC[0]+72.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiUx[1]+(90.0*phiUx[0]-210.0*phiC[0]+480.0*bcVals[0])*rdxCp2Sq[0])*omega+18.0*phiC[0]*rdxCp2Sq[1]+177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+210.0*phiC[0]*rdxCp2Sq[0])/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[1]+17.32050807568877*rdxCp2[0]*rho[0])*omega*volFac+(((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1]-36.0*phiC[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+27.0*rdxCp2[0]*phiLy[1]-354.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUy[0]+51.96152422706631*phiUx[0]+25.98076211353316*phiLy[0]-415.6921938165305*bcVals[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiUx[1]-420.0*rdxCp2Sq[0]*phiC[1]+(103.9230484541326*phiUx[0]-415.6921938165305*bcVals[0])*rdxCp2Sq[0])*omega+36.0*phiC[1]*rdxCp2Sq[1]+354.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+420.0*rdxCp2Sq[0]*phiC[1])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = ((1.732050807568877*rdxCp2[0]*rho[3]+(40.0*rdxCp2[1]+50.0*rdxCp2[0])*rho[2])*omega*volFac+((-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUy[3])+((-138.5640646055102*rdxCp2[0]*rdxCp2[1])-207.8460969082653*rdxCp2Sq[0])*phiUx[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-400.0*rdxCp2Sq[1])-500.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(120.0*rdxCp2[0]*rdxCp2[1]+180.0*rdxCp2Sq[0])*phiUx[2]+((-400.0*rdxCp2Sq[1])-500.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-1600.0*rdxCp2Sq[1])-2360.0*rdxCp2[0]*rdxCp2[1]-420.0*rdxCp2Sq[0])*phiC[2]+(346.4101615137754*phiUy[0]-346.4101615137754*phiLy[0])*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUy[1]-30.0*rdxCp2[0]*phiLy[1]+(433.0127018922193*phiUy[0]-433.0127018922193*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0])*phiC[2])/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[3] = (((40.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[3]+17.32050807568877*rdxCp2[0]*rho[2])*omega*volFac+(((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-800.0*rdxCp2[0]*rdxCp2[1])-240.0*rdxCp2Sq[0])*phiUx[3]+((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-3200.0*rdxCp2Sq[1])-4720.0*rdxCp2[0]*rdxCp2[1]-840.0*rdxCp2Sq[0])*phiC[3]-173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(692.8203230275509*rdxCp2[0]*rdxCp2[1]+207.8460969082653*rdxCp2Sq[0])*phiUx[2]-173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(692.8203230275509*phiUy[1]-692.8203230275509*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(150.0*phiUy[0]-150.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(3200.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+840.0*rdxCp2Sq[0])*phiC[3])/(3200.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+840.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((10.39230484541326*rdxCp2[0]*rho[1]-18.0*rho[0]*rdxCp2[1]-60.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-27.0*phiUy[0])-27.0*phiLy[0]+54.0*phiC[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+41.56921938165305*rdxCp2[0]*phiUx[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+((-90.0*phiUy[0])-36.0*phiUx[0]-90.0*phiLy[0]+216.0*phiC[0]+24.0*bcVals[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiUx[1]+((-60.0*phiUx[0])+60.0*phiC[0]+160.0*bcVals[0])*rdxCp2Sq[0])*omega-54.0*phiC[0]*rdxCp2Sq[1]-216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-60.0*phiC[0]*rdxCp2Sq[0]))/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[1]-17.32050807568877*rdxCp2[0]*rho[0])*omega*volFac+(((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]-60.0*rdxCp2[0]*phiUx[1]+18.0*rdxCp2[0]*phiLy[1]-216.0*rdxCp2[0]*phiC[1]+((-25.98076211353316*phiUy[0])+51.96152422706631*phiUx[0]-25.98076211353316*phiLy[0]+69.28203230275508*bcVals[0])*rdxCp2[0])*rdxCp2[1]-60.0*rdxCp2Sq[0]*phiC[1]+69.28203230275508*bcVals[0]*rdxCp2Sq[0])*omega+54.0*phiC[1]*rdxCp2Sq[1]+216.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+60.0*rdxCp2Sq[0]*phiC[1])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((5.196152422706631*rdxCp2[0]*rho[3]+((-60.0*rdxCp2[1])-30.0*rdxCp2[0])*rho[2])*omega*volFac+((-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUy[3])+(277.1281292110203*rdxCp2[0]*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0])*phiUx[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(600.0*rdxCp2Sq[1]+300.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-240.0*rdxCp2[0]*rdxCp2[1])-60.0*rdxCp2Sq[0])*phiUx[2]+(600.0*rdxCp2Sq[1]+300.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiC[2]+(519.6152422706631*phiLy[0]-519.6152422706631*phiUy[0])*rdxCp2Sq[1]+(90.0*rdxCp2[0]*phiUy[1]-90.0*rdxCp2[0]*phiLy[1]+(259.8076211353315*phiLy[0]-259.8076211353315*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-2400.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-60.0*rdxCp2Sq[0])*phiC[2]))/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[3] = (((30.0*rdxCp2[1]+3.0*rdxCp2[0])*rho[3]-8.660254037844386*rdxCp2[0]*rho[2])*omega*volFac+(((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiUx[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-60.0*rdxCp2Sq[0])*phiC[3]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]+346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[2]+86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiLy[0]-75.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiC[3])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((1.732050807568877*rdxCp2[0]*rho[1]-6.0*rho[0]*rdxCp2[1]-50.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(10.39230484541326*rdxCp2Sq[1]+86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-10.39230484541326*rdxCp2Sq[1])-86.60254037844386*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-9.0*phiUy[0])-9.0*phiLy[0]+18.0*phiC[0])*rdxCp2Sq[1]+(5.196152422706631*rdxCp2[0]*phiUy[1]+5.196152422706631*rdxCp2[0]*phiLy[1]-10.39230484541326*rdxCp2[0]*phiLx[1]-72.0*rdxCp2[0]*bcVals[1]+((-75.0*phiUy[0])-75.0*phiLy[0]-9.0*phiLx[0]+177.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-103.9230484541326*rdxCp2Sq[0]*phiLx[1]-480.0*rdxCp2Sq[0]*bcVals[1]+(210.0*phiC[0]-90.0*phiLx[0])*rdxCp2Sq[0])*omega-18.0*phiC[0]*rdxCp2Sq[1]-177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-210.0*phiC[0]*rdxCp2Sq[0]))/(18.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+210.0*rdxCp2Sq[0]); 
  phiC[1] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[1]-17.32050807568877*rdxCp2[0]*rho[0])*omega*volFac+(((-20.78460969082652*rdxCp2Sq[1])-31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(20.78460969082652*rdxCp2Sq[1]+31.17691453623978*rdxCp2[0]*rdxCp2[1])*phiLy[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*phiUy[1]+18.0*phiLy[1]-36.0*phiC[1])*rdxCp2Sq[1]+(27.0*rdxCp2[0]*phiUy[1]+27.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]-354.0*rdxCp2[0]*phiC[1]+415.6921938165305*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiUy[0])-25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]-120.0*rdxCp2Sq[0]*phiLx[1]-420.0*rdxCp2Sq[0]*phiC[1]+415.6921938165305*rdxCp2Sq[0]*bcVals[1]-103.9230484541326*phiLx[0]*rdxCp2Sq[0])*omega+36.0*phiC[1]*rdxCp2Sq[1]+354.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+420.0*rdxCp2Sq[0]*phiC[1])/(36.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[2] = -(1.0*((1.732050807568877*rdxCp2[0]*rho[3]+((-40.0*rdxCp2[1])-50.0*rdxCp2[0])*rho[2])*omega*volFac+((-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUy[3])-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-138.5640646055102*rdxCp2[0]*rdxCp2[1])-207.8460969082653*rdxCp2Sq[0])*phiLx[3]+(400.0*rdxCp2Sq[1]+500.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(400.0*rdxCp2Sq[1]+500.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+((-120.0*rdxCp2[0]*rdxCp2[1])-180.0*rdxCp2Sq[0])*phiLx[2]+(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0])*phiC[2]+(346.4101615137754*phiLy[0]-346.4101615137754*phiUy[0])*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUy[1]-30.0*rdxCp2[0]*phiLy[1]+(433.0127018922193*phiLy[0]-433.0127018922193*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+((-1600.0*rdxCp2Sq[1])-2360.0*rdxCp2[0]*rdxCp2[1]-420.0*rdxCp2Sq[0])*phiC[2]))/(1600.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+420.0*rdxCp2Sq[0]); 
  phiC[3] = (((40.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[3]-17.32050807568877*rdxCp2[0]*rho[2])*omega*volFac+(((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-800.0*rdxCp2Sq[1])-180.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-800.0*rdxCp2[0]*rdxCp2[1])-240.0*rdxCp2Sq[0])*phiLx[3]+((-3200.0*rdxCp2Sq[1])-4720.0*rdxCp2[0]*rdxCp2[1]-840.0*rdxCp2Sq[0])*phiC[3]+173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiUy[2]+173.2050807568877*rdxCp2[0]*rdxCp2[1]*phiLy[2]+((-692.8203230275509*rdxCp2[0]*rdxCp2[1])-207.8460969082653*rdxCp2Sq[0])*phiLx[2]+(692.8203230275509*phiUy[1]-692.8203230275509*phiLy[1])*rdxCp2Sq[1]+(155.8845726811989*rdxCp2[0]*phiUy[1]-155.8845726811989*rdxCp2[0]*phiLy[1]+(150.0*phiLy[0]-150.0*phiUy[0])*rdxCp2[0])*rdxCp2[1])*omega+(3200.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+840.0*rdxCp2Sq[0])*phiC[3])/(3200.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+840.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((10.39230484541326*rdxCp2[0]*rho[1]+18.0*rho[0]*rdxCp2[1]+60.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-31.17691453623978*rdxCp2Sq[1])-103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(31.17691453623978*rdxCp2Sq[1]+103.9230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*phiUy[0]+27.0*phiLy[0]-54.0*phiC[0])*rdxCp2Sq[1]+(31.17691453623978*rdxCp2[0]*phiUy[1]+31.17691453623978*rdxCp2[0]*phiLy[1]+41.56921938165305*rdxCp2[0]*phiLx[1]+24.0*rdxCp2[0]*bcVals[1]+(90.0*phiUy[0]+90.0*phiLy[0]+36.0*phiLx[0]-216.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0]*phiLx[1]+160.0*rdxCp2Sq[0]*bcVals[1]+(60.0*phiLx[0]-60.0*phiC[0])*rdxCp2Sq[0])*omega+54.0*phiC[0]*rdxCp2Sq[1]+216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+60.0*phiC[0]*rdxCp2Sq[0])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[1] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[1]+17.32050807568877*rdxCp2[0]*rho[0])*omega*volFac+(((-31.17691453623978*rdxCp2Sq[1])-20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiUy[3]+(31.17691453623978*rdxCp2Sq[1]+20.78460969082652*rdxCp2[0]*rdxCp2[1])*phiLy[3]-30.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+30.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2Sq[1]+(18.0*rdxCp2[0]*phiUy[1]+18.0*rdxCp2[0]*phiLy[1]-60.0*rdxCp2[0]*phiLx[1]-216.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiUy[0]+25.98076211353316*phiLy[0]-51.96152422706631*phiLx[0])*rdxCp2[0])*rdxCp2[1]-60.0*rdxCp2Sq[0]*phiC[1]+69.28203230275508*rdxCp2Sq[0]*bcVals[1])*omega+54.0*phiC[1]*rdxCp2Sq[1]+216.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+60.0*rdxCp2Sq[0]*phiC[1])/(54.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[2] = ((5.196152422706631*rdxCp2[0]*rho[3]+(60.0*rdxCp2[1]+30.0*rdxCp2[0])*rho[2])*omega*volFac+((-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUy[3])-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLy[3]+(277.1281292110203*rdxCp2[0]*rdxCp2[1]+69.28203230275508*rdxCp2Sq[0])*phiLx[3]+((-600.0*rdxCp2Sq[1])-300.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+((-600.0*rdxCp2Sq[1])-300.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(240.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiLx[2]+((-2400.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-60.0*rdxCp2Sq[0])*phiC[2]+(519.6152422706631*phiUy[0]-519.6152422706631*phiLy[0])*rdxCp2Sq[1]+(90.0*rdxCp2[0]*phiUy[1]-90.0*rdxCp2[0]*phiLy[1]+(259.8076211353315*phiUy[0]-259.8076211353315*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiC[2])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 
  phiC[3] = (((30.0*rdxCp2[1]+3.0*rdxCp2[0])*rho[3]+8.660254037844386*rdxCp2[0]*rho[2])*omega*volFac+(((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-2400.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-60.0*rdxCp2Sq[0])*phiC[3]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiUy[2]-86.60254037844386*rdxCp2[0]*rdxCp2[1]*phiLy[2]-346.4101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(519.6152422706631*phiUy[1]-519.6152422706631*phiLy[1])*rdxCp2Sq[1]+(51.96152422706631*rdxCp2[0]*phiUy[1]-51.96152422706631*rdxCp2[0]*phiLy[1]+(75.0*phiUy[0]-75.0*phiLy[0])*rdxCp2[0])*rdxCp2[1])*omega+(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0])*phiC[3])/(2400.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+60.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((1.732050807568877*rdxCp2[1]*rho[2]+50.0*rho[0]*rdxCp2[1]+6.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiUy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(480.0*rdxCp2Sq[1]+72.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(90.0*phiUy[0]-210.0*phiC[0])*rdxCp2Sq[1]+((-86.60254037844386*rdxCp2[0]*phiUx[1])+86.60254037844386*rdxCp2[0]*phiLx[1]+(9.0*phiUy[0]+75.0*phiUx[0]+75.0*phiLx[0]-177.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-10.39230484541326*rdxCp2Sq[0]*phiUx[1]+10.39230484541326*rdxCp2Sq[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2Sq[0])*omega+210.0*phiC[0]*rdxCp2Sq[1]+177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+18.0*phiC[0]*rdxCp2Sq[0])/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = ((1.732050807568877*rdxCp2[1]*rho[3]+(50.0*rdxCp2[1]+40.0*rdxCp2[0])*rho[1])*omega*volFac+(((-207.8460969082653*rdxCp2Sq[1])-138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiUy[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(180.0*phiUy[1]-420.0*phiC[1])*rdxCp2Sq[1]+(120.0*rdxCp2[0]*phiUy[1]-500.0*rdxCp2[0]*phiUx[1]-500.0*rdxCp2[0]*phiLx[1]-2360.0*rdxCp2[0]*phiC[1]+(433.0127018922193*phiUx[0]-433.0127018922193*phiLx[0])*rdxCp2[0])*rdxCp2[1]-400.0*rdxCp2Sq[0]*phiUx[1]-400.0*rdxCp2Sq[0]*phiLx[1]-1600.0*rdxCp2Sq[0]*phiC[1]+(346.4101615137754*phiUx[0]-346.4101615137754*phiLx[0])*rdxCp2Sq[0])*omega+420.0*phiC[1]*rdxCp2Sq[1]+2360.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+1600.0*rdxCp2Sq[0]*phiC[1])/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 
  phiC[2] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[2]+17.32050807568877*rho[0]*rdxCp2[1])*omega*volFac+(((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiUy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]+((-420.0*rdxCp2Sq[1])-354.0*rdxCp2[0]*rdxCp2[1]-36.0*rdxCp2Sq[0])*phiC[2]+((-415.6921938165305*rdxCp2Sq[1])-415.6921938165305*rdxCp2[0]*rdxCp2[1])*bcVals[2]+103.9230484541326*phiUy[0]*rdxCp2Sq[1]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]+25.98076211353316*phiUx[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0])*phiC[2])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((9.0*rdxCp2[1]+40.0*rdxCp2[0])*rho[3]+17.32050807568877*rdxCp2[1]*rho[1])*omega*volFac+(((-240.0*rdxCp2Sq[1])-800.0*rdxCp2[0]*rdxCp2[1])*phiUy[3]+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiUx[3]+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiLx[3]+((-840.0*rdxCp2Sq[1])-4720.0*rdxCp2[0]*rdxCp2[1]-3200.0*rdxCp2Sq[0])*phiC[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-692.8203230275509*rdxCp2Sq[0])*phiLx[2]+207.8460969082653*phiUy[1]*rdxCp2Sq[1]+(692.8203230275509*rdxCp2[0]*phiUy[1]-173.2050807568877*rdxCp2[0]*phiUx[1]-173.2050807568877*rdxCp2[0]*phiLx[1]+(150.0*phiUx[0]-150.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(840.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+3200.0*rdxCp2Sq[0])*phiC[3])/(840.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+3200.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((10.39230484541326*rdxCp2[1]*rho[2]-60.0*rho[0]*rdxCp2[1]-18.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiUy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(60.0*phiC[0]-60.0*phiUy[0])*rdxCp2Sq[1]+(103.9230484541326*rdxCp2[0]*phiUx[1]-103.9230484541326*rdxCp2[0]*phiLx[1]+((-36.0*phiUy[0])-90.0*phiUx[0]-90.0*phiLx[0]+216.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0]*phiUx[1]-31.17691453623978*rdxCp2Sq[0]*phiLx[1]+((-27.0*phiUx[0])-27.0*phiLx[0]+54.0*phiC[0])*rdxCp2Sq[0])*omega-60.0*phiC[0]*rdxCp2Sq[1]-216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-54.0*phiC[0]*rdxCp2Sq[0]))/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((5.196152422706631*rdxCp2[1]*rho[3]+((-30.0*rdxCp2[1])-60.0*rdxCp2[0])*rho[1])*omega*volFac+((69.28203230275508*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*phiUy[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUx[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLx[3]+90.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-90.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(60.0*phiC[1]-60.0*phiUy[1])*rdxCp2Sq[1]+((-240.0*rdxCp2[0]*phiUy[1])+300.0*rdxCp2[0]*phiUx[1]+300.0*rdxCp2[0]*phiLx[1]+1440.0*rdxCp2[0]*phiC[1]+(259.8076211353315*phiLx[0]-259.8076211353315*phiUx[0])*rdxCp2[0])*rdxCp2[1]+600.0*rdxCp2Sq[0]*phiUx[1]+600.0*rdxCp2Sq[0]*phiLx[1]+2400.0*rdxCp2Sq[0]*phiC[1]+(519.6152422706631*phiLx[0]-519.6152422706631*phiUx[0])*rdxCp2Sq[0])*omega-60.0*phiC[1]*rdxCp2Sq[1]-1440.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-2400.0*rdxCp2Sq[0]*phiC[1]))/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 
  phiC[2] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[2]-17.32050807568877*rho[0]*rdxCp2[1])*omega*volFac+(((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]-60.0*rdxCp2[0]*rdxCp2[1]*phiUy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+((-60.0*rdxCp2Sq[1])-216.0*rdxCp2[0]*rdxCp2[1]-54.0*rdxCp2Sq[0])*phiC[2]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[2]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+(51.96152422706631*phiUy[0]-25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0])*phiC[2])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((3.0*rdxCp2[1]+30.0*rdxCp2[0])*rho[3]-8.660254037844386*rdxCp2[1]*rho[1])*omega*volFac+((-400.0*rdxCp2[0]*rdxCp2[1]*phiUy[3])+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+((-60.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-2400.0*rdxCp2Sq[0])*phiC[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+(346.4101615137754*rdxCp2[0]*phiUy[1]+86.60254037844386*rdxCp2[0]*phiUx[1]+86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiLx[0]-75.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0])*phiC[3])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((1.732050807568877*rdxCp2[1]*rho[2]-50.0*rho[0]*rdxCp2[1]-6.0*rdxCp2[0]*rho[0])*omega*volFac+((-6.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+6.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+((-480.0*rdxCp2Sq[1])-72.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiUx[2]+((-103.9230484541326*rdxCp2Sq[1])-10.39230484541326*rdxCp2[0]*rdxCp2[1])*phiLy[2]+5.196152422706631*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(210.0*phiC[0]-90.0*phiLy[0])*rdxCp2Sq[1]+(86.60254037844386*rdxCp2[0]*phiUx[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+((-75.0*phiUx[0])-9.0*phiLy[0]-75.0*phiLx[0]+177.0*phiC[0])*rdxCp2[0])*rdxCp2[1]+10.39230484541326*rdxCp2Sq[0]*phiUx[1]-10.39230484541326*rdxCp2Sq[0]*phiLx[1]+((-9.0*phiUx[0])-9.0*phiLx[0]+18.0*phiC[0])*rdxCp2Sq[0])*omega-210.0*phiC[0]*rdxCp2Sq[1]-177.0*phiC[0]*rdxCp2[0]*rdxCp2[1]-18.0*phiC[0]*rdxCp2Sq[0]))/(210.0*rdxCp2Sq[1]+177.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0]); 
  phiC[1] = -(1.0*((1.732050807568877*rdxCp2[1]*rho[3]+((-50.0*rdxCp2[1])-40.0*rdxCp2[0])*rho[1])*omega*volFac+((-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiUx[3])+((-207.8460969082653*rdxCp2Sq[1])-138.5640646055102*rdxCp2[0]*rdxCp2[1])*phiLy[3]-34.64101615137754*rdxCp2[0]*rdxCp2[1]*phiLx[3]+30.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-30.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(420.0*phiC[1]-180.0*phiLy[1])*rdxCp2Sq[1]+(500.0*rdxCp2[0]*phiUx[1]-120.0*rdxCp2[0]*phiLy[1]+500.0*rdxCp2[0]*phiLx[1]+2360.0*rdxCp2[0]*phiC[1]+(433.0127018922193*phiLx[0]-433.0127018922193*phiUx[0])*rdxCp2[0])*rdxCp2[1]+400.0*rdxCp2Sq[0]*phiUx[1]+400.0*rdxCp2Sq[0]*phiLx[1]+1600.0*rdxCp2Sq[0]*phiC[1]+(346.4101615137754*phiLx[0]-346.4101615137754*phiUx[0])*rdxCp2Sq[0])*omega-420.0*phiC[1]*rdxCp2Sq[1]-2360.0*rdxCp2[0]*phiC[1]*rdxCp2[1]-1600.0*rdxCp2Sq[0]*phiC[1]))/(420.0*rdxCp2Sq[1]+2360.0*rdxCp2[0]*rdxCp2[1]+1600.0*rdxCp2Sq[0]); 
  phiC[2] = (((9.0*rdxCp2[1]+6.0*rdxCp2[0])*rho[2]-17.32050807568877*rho[0]*rdxCp2[1])*omega*volFac+(((-31.17691453623978*rdxCp2[0]*rdxCp2[1])-20.78460969082652*rdxCp2Sq[0])*phiUx[3]+(31.17691453623978*rdxCp2[0]*rdxCp2[1]+20.78460969082652*rdxCp2Sq[0])*phiLx[3]+(415.6921938165305*rdxCp2Sq[1]+415.6921938165305*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiUx[2]+((-120.0*rdxCp2Sq[1])-60.0*rdxCp2[0]*rdxCp2[1])*phiLy[2]+(27.0*rdxCp2[0]*rdxCp2[1]+18.0*rdxCp2Sq[0])*phiLx[2]+((-420.0*rdxCp2Sq[1])-354.0*rdxCp2[0]*rdxCp2[1]-36.0*rdxCp2Sq[0])*phiC[2]-103.9230484541326*phiLy[0]*rdxCp2Sq[1]+(30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]+((-25.98076211353316*phiUx[0])-51.96152422706631*phiLy[0]-25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0])*phiC[2])/(420.0*rdxCp2Sq[1]+354.0*rdxCp2[0]*rdxCp2[1]+36.0*rdxCp2Sq[0]); 
  phiC[3] = (((9.0*rdxCp2[1]+40.0*rdxCp2[0])*rho[3]-17.32050807568877*rdxCp2[1]*rho[1])*omega*volFac+(((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiUx[3]+((-240.0*rdxCp2Sq[1])-800.0*rdxCp2[0]*rdxCp2[1])*phiLy[3]+((-180.0*rdxCp2[0]*rdxCp2[1])-800.0*rdxCp2Sq[0])*phiLx[3]+((-840.0*rdxCp2Sq[1])-4720.0*rdxCp2[0]*rdxCp2[1]-3200.0*rdxCp2Sq[0])*phiC[3]+(155.8845726811989*rdxCp2[0]*rdxCp2[1]+692.8203230275509*rdxCp2Sq[0])*phiUx[2]+((-155.8845726811989*rdxCp2[0]*rdxCp2[1])-692.8203230275509*rdxCp2Sq[0])*phiLx[2]-207.8460969082653*phiLy[1]*rdxCp2Sq[1]+(173.2050807568877*rdxCp2[0]*phiUx[1]-692.8203230275509*rdxCp2[0]*phiLy[1]+173.2050807568877*rdxCp2[0]*phiLx[1]+(150.0*phiLx[0]-150.0*phiUx[0])*rdxCp2[0])*rdxCp2[1])*omega+(840.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+3200.0*rdxCp2Sq[0])*phiC[3])/(840.0*rdxCp2Sq[1]+4720.0*rdxCp2[0]*rdxCp2[1]+3200.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((10.39230484541326*rdxCp2[1]*rho[2]+60.0*rho[0]*rdxCp2[1]+18.0*rdxCp2[0]*rho[0])*omega*volFac+((-36.0*rdxCp2[0]*rdxCp2[1]*phiUx[3])+36.0*rdxCp2[0]*rdxCp2[1]*phiLx[3]+(160.0*rdxCp2Sq[1]+24.0*rdxCp2[0]*rdxCp2[1])*bcVals[3]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiUx[2]+(69.28203230275508*rdxCp2Sq[1]+41.56921938165305*rdxCp2[0]*rdxCp2[1])*phiLy[2]+31.17691453623978*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(60.0*phiLy[0]-60.0*phiC[0])*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2[0]*phiUx[1])+103.9230484541326*rdxCp2[0]*phiLx[1]+(90.0*phiUx[0]+36.0*phiLy[0]+90.0*phiLx[0]-216.0*phiC[0])*rdxCp2[0])*rdxCp2[1]-31.17691453623978*rdxCp2Sq[0]*phiUx[1]+31.17691453623978*rdxCp2Sq[0]*phiLx[1]+(27.0*phiUx[0]+27.0*phiLx[0]-54.0*phiC[0])*rdxCp2Sq[0])*omega+60.0*phiC[0]*rdxCp2Sq[1]+216.0*phiC[0]*rdxCp2[0]*rdxCp2[1]+54.0*phiC[0]*rdxCp2Sq[0])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[1] = ((5.196152422706631*rdxCp2[1]*rho[3]+(30.0*rdxCp2[1]+60.0*rdxCp2[0])*rho[1])*omega*volFac+((-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiUx[3])+(69.28203230275508*rdxCp2Sq[1]+277.1281292110203*rdxCp2[0]*rdxCp2[1])*phiLy[3]-103.9230484541326*rdxCp2[0]*rdxCp2[1]*phiLx[3]+90.0*rdxCp2[0]*rdxCp2[1]*phiUx[2]-90.0*rdxCp2[0]*rdxCp2[1]*phiLx[2]+(60.0*phiLy[1]-60.0*phiC[1])*rdxCp2Sq[1]+((-300.0*rdxCp2[0]*phiUx[1])+240.0*rdxCp2[0]*phiLy[1]-300.0*rdxCp2[0]*phiLx[1]-1440.0*rdxCp2[0]*phiC[1]+(259.8076211353315*phiUx[0]-259.8076211353315*phiLx[0])*rdxCp2[0])*rdxCp2[1]-600.0*rdxCp2Sq[0]*phiUx[1]-600.0*rdxCp2Sq[0]*phiLx[1]-2400.0*rdxCp2Sq[0]*phiC[1]+(519.6152422706631*phiUx[0]-519.6152422706631*phiLx[0])*rdxCp2Sq[0])*omega+60.0*phiC[1]*rdxCp2Sq[1]+1440.0*rdxCp2[0]*phiC[1]*rdxCp2[1]+2400.0*rdxCp2Sq[0]*phiC[1])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 
  phiC[2] = (((6.0*rdxCp2[1]+9.0*rdxCp2[0])*rho[2]+17.32050807568877*rho[0]*rdxCp2[1])*omega*volFac+(((-20.78460969082652*rdxCp2[0]*rdxCp2[1])-31.17691453623978*rdxCp2Sq[0])*phiUx[3]+(20.78460969082652*rdxCp2[0]*rdxCp2[1]+31.17691453623978*rdxCp2Sq[0])*phiLx[3]+(69.28203230275508*rdxCp2Sq[1]+69.28203230275508*rdxCp2[0]*rdxCp2[1])*bcVals[3]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiUx[2]-60.0*rdxCp2[0]*rdxCp2[1]*phiLy[2]+(18.0*rdxCp2[0]*rdxCp2[1]+27.0*rdxCp2Sq[0])*phiLx[2]+((-60.0*rdxCp2Sq[1])-216.0*rdxCp2[0]*rdxCp2[1]-54.0*rdxCp2Sq[0])*phiC[2]+((-30.0*rdxCp2[0]*phiUx[1])+30.0*rdxCp2[0]*phiLx[1]+(25.98076211353316*phiUx[0]-51.96152422706631*phiLy[0]+25.98076211353316*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0])*phiC[2])/(60.0*rdxCp2Sq[1]+216.0*rdxCp2[0]*rdxCp2[1]+54.0*rdxCp2Sq[0]); 
  phiC[3] = (((3.0*rdxCp2[1]+30.0*rdxCp2[0])*rho[3]+8.660254037844386*rdxCp2[1]*rho[1])*omega*volFac+(((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiUx[3]-400.0*rdxCp2[0]*rdxCp2[1]*phiLy[3]+((-60.0*rdxCp2[0]*rdxCp2[1])-600.0*rdxCp2Sq[0])*phiLx[3]+((-60.0*rdxCp2Sq[1])-1440.0*rdxCp2[0]*rdxCp2[1]-2400.0*rdxCp2Sq[0])*phiC[3]+(51.96152422706631*rdxCp2[0]*rdxCp2[1]+519.6152422706631*rdxCp2Sq[0])*phiUx[2]+((-51.96152422706631*rdxCp2[0]*rdxCp2[1])-519.6152422706631*rdxCp2Sq[0])*phiLx[2]+((-86.60254037844386*rdxCp2[0]*phiUx[1])-346.4101615137754*rdxCp2[0]*phiLy[1]-86.60254037844386*rdxCp2[0]*phiLx[1]+(75.0*phiUx[0]-75.0*phiLx[0])*rdxCp2[0])*rdxCp2[1])*omega+(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0])*phiC[3])/(60.0*rdxCp2Sq[1]+1440.0*rdxCp2[0]*rdxCp2[1]+2400.0*rdxCp2Sq[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(727.4613391789284*rdxCp2Sq[1]+4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+(4382.08854314926*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-43647.6803507357*rdxCp2R3[1])-269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(2182.384017536785*rdxCp2[0]*rdxCp2Sq[1]+16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(201600.0*rdxCp2R3[1]+1259760.0*rdxCp2[0]*rdxCp2Sq[1]+252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(37800.0*phiUy[0]-88200.0*phiC[0])*rdxCp2R3[1]+(16211.99555884469*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiUx[1]+(233370.0*phiUy[0]+31500.0*phiUx[0]-642810.0*phiC[0]+252000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(2182.384017536785*rdxCp2Sq[0]*phiUy[1]-269472.4646415659*rdxCp2Sq[0]*phiUx[1]+(31500.0*phiUy[0]+233370.0*phiUx[0]-642810.0*phiC[0]+1259760.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-43647.6803507357*rdxCp2R3[0]*phiUx[1]+(37800.0*phiUx[0]-88200.0*phiC[0]+201600.0*bcVals[0])*rdxCp2R3[0])*omega+88200.0*phiC[0]*rdxCp2R3[1]+642810.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+642810.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+88200.0*phiC[0]*rdxCp2R3[0])/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = (((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]+7274.613391789284*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-93600.0*rdxCp2[0]*rdxCp2Sq[1])-12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(12600.0*rdxCp2[0]*rdxCp2Sq[1]+8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(403221.4280020346*rdxCp2[0]*rdxCp2Sq[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(75600.0*phiUy[1]-176400.0*phiC[1])*rdxCp2R3[1]+(82980.0*rdxCp2[0]*phiUy[1]-210000.0*rdxCp2[0]*phiUx[1]-1285620.0*rdxCp2[0]*phiC[1]+(81059.97779422344*phiUy[0]+181865.3347947321*phiUx[0]-1454922.678357857*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(11340.0*rdxCp2Sq[0]*phiUy[1]-341400.0*rdxCp2Sq[0]*phiUx[1]-1285620.0*rdxCp2Sq[0]*phiC[1]+(10911.92008768392*phiUy[0]+295661.0728520073*phiUx[0]-1313587.332460236*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-50400.0*rdxCp2R3[0]*phiUx[1]-176400.0*rdxCp2R3[0]*phiC[1]+(43647.6803507357*phiUx[0]-174590.7214029428*bcVals[0])*rdxCp2R3[0])*omega+176400.0*phiC[1]*rdxCp2R3[1]+1285620.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+176400.0*rdxCp2R3[0]*phiC[1])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = (((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+7274.613391789284*rho[0]*rdxCp2Sq[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiUx[3]+((-50400.0*rdxCp2R3[1])-341400.0*rdxCp2[0]*rdxCp2Sq[1]-210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(11340.0*rdxCp2[0]*rdxCp2Sq[1]+82980.0*rdxCp2Sq[0]*rdxCp2[1]+75600.0*rdxCp2R3[0])*phiUx[2]+((-176400.0*rdxCp2R3[1])-1285620.0*rdxCp2[0]*rdxCp2Sq[1]-1285620.0*rdxCp2Sq[0]*rdxCp2[1]-176400.0*rdxCp2R3[0])*phiC[2]+((-174590.7214029428*rdxCp2R3[1])-1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]-1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+43647.6803507357*phiUy[0]*rdxCp2R3[1]+(8640.0*rdxCp2[0]*phiUy[1]-12600.0*rdxCp2[0]*phiUx[1]+(295661.0728520073*phiUy[0]+10911.92008768392*phiUx[0]+87295.36070147139*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(12600.0*rdxCp2Sq[0]*phiUy[1]-93600.0*rdxCp2Sq[0]*phiUx[1]+(181865.3347947321*phiUy[0]+81059.97779422344*phiUx[0]+403221.4280020346*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0])*phiC[2])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(640.8587988004846*rdxCp2[0]*rdxCp2[1]+2424.871130596428*rdxCp2Sq[0])*rho[2]+(2424.871130596428*rdxCp2Sq[1]+640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]+5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiUx[3]+((-117600.0*rdxCp2R3[1])-857080.0*rdxCp2[0]*rdxCp2Sq[1]-857080.0*rdxCp2Sq[0]*rdxCp2[1]-117600.0*rdxCp2R3[0])*phiC[3]+((-16627.68775266122*rdxCp2[0]*rdxCp2Sq[1])-24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(21823.84017536785*rdxCp2[0]*rdxCp2Sq[1]+128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]+29098.45356715714*rdxCp2R3[0])*phiUx[2]+(26400.0*rdxCp2[0]*rdxCp2Sq[1]-168000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+29098.45356715714*phiUy[1]*rdxCp2R3[1]+(128933.8621154272*rdxCp2[0]*phiUy[1]-24248.71130596428*rdxCp2[0]*phiUx[1]+(14400.0*phiUy[0]+21000.0*phiUx[0]-168000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]-16627.68775266122*rdxCp2Sq[0]*phiUx[1]+(21000.0*phiUy[0]+14400.0*phiUx[0]+26400.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0])*phiC[3])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(207.8460969082653*rdxCp2R3[1]+6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]+13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1143.153532995459*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-1091.192008768392*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+((2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(1385.640646055102*rdxCp2R4[1]+42816.29596310264*rdxCp2[0]*rdxCp2R3[1]+122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]+72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(623.5382907247956*rdxCp2[0]*rdxCp2R3[1]+22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]+48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(3200.0*rdxCp2R4[1]+96720.0*rdxCp2[0]*rdxCp2R3[1]+222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(1200.0*phiC[0]-1200.0*phiUy[0])*rdxCp2R4[1]+((-2078.460969082652*rdxCp2[0]*phiUy[1])+2078.460969082652*rdxCp2[0]*phiUx[1]+((-37080.0*phiUy[0])-1800.0*phiUx[0]+42480.0*phiC[0]-14400.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-6131.459858793824*rdxCp2Sq[0]*phiUy[1])+74720.67183852136*rdxCp2Sq[0]*phiUx[1]+((-106140.0*phiUy[0])-64710.0*phiUx[0]+260670.0*phiC[0]-359280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-4364.768035073569*rdxCp2R3[0]*phiUy[1])+188308.5637988883*rdxCp2R3[0]*phiUx[1]+((-63000.0*phiUy[0])-163080.0*phiUx[0]+446040.0*phiC[0]-879840.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+65471.52052610354*rdxCp2R4[0]*phiUx[1]+((-56700.0*phiUx[0])+132300.0*phiC[0]-302400.0*bcVals[0])*rdxCp2R4[0])*omega-1200.0*phiC[0]*rdxCp2R4[1]-42480.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-446040.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-132300.0*phiC[0]*rdxCp2R4[0]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]-11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-10911.92008768392*rdxCp2R3[0]*rho[0])*omega*volFac+((2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(12000.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25200.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(3600.0*rdxCp2[0]*rdxCp2R3[1]+25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(31869.73485926734*rdxCp2[0]*rdxCp2R3[1]+81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]+14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(2400.0*phiC[1]-2400.0*phiUy[1])*rdxCp2R4[1]+((-24960.0*rdxCp2[0]*phiUy[1])+12000.0*rdxCp2[0]*phiUx[1]+84960.0*rdxCp2[0]*phiC[1]+((-10392.30484541326*phiUy[0])-10392.30484541326*phiUx[0]+83138.4387633061*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-67140.0*rdxCp2Sq[0]*phiUy[1])+114600.0*rdxCp2Sq[0]*phiUx[1]+521340.0*rdxCp2Sq[0]*phiC[1]+((-30657.29929396912*phiUy[0])-99246.51127369665*phiUx[0]+519615.2422706631*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-22680.0*rdxCp2R3[0]*phiUy[1])+237600.0*rdxCp2R3[0]*phiUx[1]+892080.0*rdxCp2R3[0]*phiC[1]+((-21823.84017536785*phiUy[0])-205767.6359391825*phiUx[0]+910365.9044582016*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+75600.0*rdxCp2R4[0]*phiUx[1]+264600.0*rdxCp2R4[0]*phiC[1]+(261886.0821044141*bcVals[0]-65471.52052610354*phiUx[0])*rdxCp2R4[0])*omega-2400.0*phiC[1]*rdxCp2R4[1]-84960.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-892080.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-264600.0*rdxCp2R4[0]*phiC[1]))/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = (((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+(120.0*rdxCp2R3[1]+3870.0*rdxCp2[0]*rdxCp2Sq[1]+15150.0*rdxCp2Sq[0]*rdxCp2[1]+15750.0*rdxCp2R3[0])*rho[2]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiUx[3]+((-1800.0*rdxCp2[0]*rdxCp2R3[1])-50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-105000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(360.0*rdxCp2[0]*rdxCp2R3[1]+12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]+50760.0*rdxCp2R3[0]*rdxCp2[1]+56700.0*rdxCp2R4[0])*phiUx[2]+((-1200.0*rdxCp2R4[1])-42480.0*rdxCp2[0]*rdxCp2R3[1]-260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]-446040.0*rdxCp2R3[0]*rdxCp2[1]-132300.0*rdxCp2R4[0])*phiC[2]+(1385.640646055102*rdxCp2R4[1]+43647.6803507357*rdxCp2[0]*rdxCp2R3[1]+145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]+121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-600.0*rdxCp2[0]*phiUy[1])+600.0*rdxCp2[0]*phiUx[1]+(1558.845726811989*phiUy[0]-519.6152422706631*phiUx[0]-4156.921938165305*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(21600.0*rdxCp2Sq[0]*phiUx[1]+(43647.6803507357*phiUy[0]-18706.14872174387*phiUx[0]-99766.1265159673*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(6300.0*rdxCp2R3[0]*phiUy[1]+46800.0*rdxCp2R3[0]*phiUx[1]+(90932.66739736605*phiUy[0]-40529.98889711172*phiUx[0]-201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0])*phiC[2])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiUx[3]+((-2400.0*rdxCp2R4[1])-84960.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892080.0*rdxCp2R3[0]*rdxCp2[1]-264600.0*rdxCp2R4[0])*phiC[3]+(3464.101615137754*rdxCp2[0]*rdxCp2R3[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(2078.460969082652*rdxCp2[0]*rdxCp2R3[1]+39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]+145907.9600296021*rdxCp2R3[0]*rdxCp2[1]+65471.52052610354*rdxCp2R4[0])*phiUx[2]+(10400.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(17320.50807568877*rdxCp2[0]*phiUy[1]+3464.101615137754*rdxCp2[0]*phiUx[1]+((-3000.0*phiUy[0])-3000.0*phiUx[0]+24000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(87295.36070147139*rdxCp2Sq[0]*phiUy[1]+24941.53162899183*rdxCp2Sq[0]*phiUx[1]+(86400.0*bcVals[0]-21600.0*phiUx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(32735.76026305177*rdxCp2R3[0]*phiUy[1]+24941.53162899183*rdxCp2R3[0]*phiUx[1]+(31500.0*phiUy[0]-21600.0*phiUx[0]-39600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0])*phiC[3])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1091.192008768392*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(13146.26562944778*rdxCp2[0]*rdxCp2Sq[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(65471.52052610354*rdxCp2R4[1]+188308.5637988883*rdxCp2[0]*rdxCp2R3[1]+74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-4364.768035073569*rdxCp2[0]*rdxCp2R3[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-302400.0*rdxCp2R4[1])-879840.0*rdxCp2[0]*rdxCp2R3[1]-359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(132300.0*phiC[0]-56700.0*phiUy[0])*rdxCp2R4[1]+(48635.98667653406*rdxCp2[0]*phiUy[1]+72746.13391789283*rdxCp2[0]*phiUx[1]+((-163080.0*phiUy[0])-63000.0*phiUx[0]+446040.0*phiC[0]+42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(22447.37846609264*rdxCp2Sq[0]*phiUy[1]+122559.9151435737*rdxCp2Sq[0]*phiUx[1]+((-64710.0*phiUy[0])-106140.0*phiUx[0]+260670.0*phiC[0]+222560.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(623.5382907247956*rdxCp2R3[0]*phiUy[1]+42816.29596310264*rdxCp2R3[0]*phiUx[1]+((-1800.0*phiUy[0])-37080.0*phiUx[0]+42480.0*phiC[0]+96720.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+1385.640646055102*rdxCp2R4[0]*phiUx[1]+((-1200.0*phiUx[0])+1200.0*phiC[0]+3200.0*bcVals[0])*rdxCp2R4[0])*omega-132300.0*phiC[0]*rdxCp2R4[1]-446040.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-42480.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-1200.0*phiC[0]*rdxCp2R4[0]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = (((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-885.0*rdxCp2[0]*rdxCp2Sq[1])-360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(15750.0*rdxCp2R3[1]+15150.0*rdxCp2[0]*rdxCp2Sq[1]+3870.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+(((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+(46800.0*rdxCp2[0]*rdxCp2R3[1]+21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(6300.0*rdxCp2[0]*rdxCp2R3[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(56700.0*phiUy[1]-132300.0*phiC[1])*rdxCp2R4[1]+(50760.0*rdxCp2[0]*phiUy[1]-105000.0*rdxCp2[0]*phiUx[1]-446040.0*rdxCp2[0]*phiC[1]+((-40529.98889711172*phiUy[0])+90932.66739736605*phiUx[0]+121243.5565298214*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(12870.0*rdxCp2Sq[0]*phiUy[1]-50400.0*rdxCp2Sq[0]*phiUx[1]-260670.0*rdxCp2Sq[0]*phiC[1]+((-18706.14872174387*phiUy[0])+43647.6803507357*phiUx[0]+145838.6779972995*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(360.0*rdxCp2R3[0]*phiUy[1]-1800.0*rdxCp2R3[0]*phiUx[1]-42480.0*rdxCp2R3[0]*phiC[1]+((-519.6152422706631*phiUy[0])+1558.845726811989*phiUx[0]+43647.6803507357*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-1200.0*rdxCp2R4[0]*phiC[1]+1385.640646055102*bcVals[0]*rdxCp2R4[0])*omega+132300.0*phiC[1]*rdxCp2R4[1]+446040.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+1200.0*rdxCp2R4[0]*phiC[1])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-10911.92008768392*rho[0]*rdxCp2R3[1]-30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiUx[3]+(75600.0*rdxCp2R4[1]+237600.0*rdxCp2[0]*rdxCp2R3[1]+114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+12000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-22680.0*rdxCp2[0]*rdxCp2R3[1])-67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiUx[2]+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[2]+(261886.0821044141*rdxCp2R4[1]+910365.9044582016*rdxCp2[0]*rdxCp2R3[1]+519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]+83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]-65471.52052610354*phiUy[0]*rdxCp2R4[1]+(25920.0*rdxCp2[0]*phiUy[1]+25200.0*rdxCp2[0]*phiUx[1]+((-205767.6359391825*phiUy[0])-21823.84017536785*phiUx[0]+14549.22678357857*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(25920.0*rdxCp2Sq[0]*phiUy[1]+35400.0*rdxCp2Sq[0]*phiUx[1]+((-99246.51127369665*phiUy[0])-30657.29929396912*phiUx[0]+81752.798117251*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3600.0*rdxCp2R3[0]*phiUy[1]+12000.0*rdxCp2R3[0]*phiUx[1]+((-10392.30484541326*phiUy[0])-10392.30484541326*phiUx[0]+31869.73485926734*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[2]))/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[3]+(24941.53162899183*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(32735.76026305177*rdxCp2[0]*rdxCp2R3[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]+17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+((-39600.0*rdxCp2[0]*rdxCp2R3[1])+86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+65471.52052610354*phiUy[1]*rdxCp2R4[1]+(145907.9600296021*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiUx[1]+((-21600.0*phiUy[0])+31500.0*phiUx[0]+42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(39386.83536411626*rdxCp2Sq[0]*phiUy[1]+(35400.0*bcVals[0]-21600.0*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(2078.460969082652*rdxCp2R3[0]*phiUy[1]+3464.101615137754*rdxCp2R3[0]*phiUx[1]+((-3000.0*phiUy[0])-3000.0*phiUx[0]+10400.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[3])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-25.98076211353316*rdxCp2Sq[1])-285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+((-285.7883832488647*rdxCp2[0]*rdxCp2[1])-25.98076211353316*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-173.2050807568877*rdxCp2R3[1])-1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-103.9230484541326*rdxCp2[0]*rdxCp2Sq[1])-519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-400.0*rdxCp2R3[1])-4440.0*rdxCp2[0]*rdxCp2Sq[1]-200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(150.0*phiUy[0]-150.0*phiC[0])*rdxCp2R3[1]+((-519.6152422706631*rdxCp2[0]*phiUy[1])-346.4101615137754*rdxCp2[0]*phiUx[1]+(1710.0*phiUy[0]+300.0*phiUx[0]-2010.0*phiC[0]-200.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2Sq[0]*phiUy[1])-1974.53792062852*rdxCp2Sq[0]*phiUx[1]+(300.0*phiUy[0]+1710.0*phiUx[0]-2010.0*phiC[0]-4440.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-173.2050807568877*rdxCp2R3[0]*phiUx[1]+(150.0*phiUx[0]-150.0*phiC[0]-400.0*bcVals[0])*rdxCp2R3[0])*omega+150.0*phiC[0]*rdxCp2R3[1]+2010.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+2010.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+150.0*phiC[0]*rdxCp2R3[0])/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0]*rho[0])*omega*volFac+((1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(900.0*rdxCp2[0]*rdxCp2Sq[1]-900.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-7967.433714816835*rdxCp2[0]*rdxCp2Sq[1])-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(900.0*phiC[1]-900.0*phiUy[1])*rdxCp2R3[1]+((-3060.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiUx[1]+12060.0*rdxCp2[0]*phiC[1]+(2598.076211353316*phiUy[0]-2598.076211353316*phiUx[0]-3464.101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-360.0*rdxCp2Sq[0]*phiUy[1])+600.0*rdxCp2Sq[0]*phiUx[1]+12060.0*rdxCp2Sq[0]*phiC[1]+(519.6152422706631*phiUy[0]-519.6152422706631*phiUx[0]-12124.35565298214*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+900.0*rdxCp2R3[0]*phiC[1]-1039.230484541326*bcVals[0]*rdxCp2R3[0])*omega-900.0*phiC[1]*rdxCp2R3[1]-12060.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-900.0*rdxCp2R3[0]*phiC[1]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+259.8076211353315*rho[0]*rdxCp2Sq[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiUx[3]+(600.0*rdxCp2[0]*rdxCp2Sq[1]+3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-3060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiUx[2]+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[2]+((-1039.230484541326*rdxCp2R3[1])-12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-900.0*rdxCp2[0]*phiUy[1])-600.0*rdxCp2[0]*phiUx[1]+((-519.6152422706631*phiUy[0])+519.6152422706631*phiUx[0]-346.4101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(900.0*rdxCp2Sq[0]*phiUy[1]-3000.0*rdxCp2Sq[0]*phiUx[1]+((-2598.076211353316*phiUy[0])+2598.076211353316*phiUx[0]-7967.433714816835*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[2]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-129.9038105676658*rdxCp2Sq[0])*rho[2]+((-129.9038105676658*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[3]+(866.0254037844386*rdxCp2Sq[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2Sq[1])*phiUy[2]+(519.6152422706631*rdxCp2[0]*rdxCp2Sq[1]+2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-2600.0*rdxCp2[0]*rdxCp2Sq[1])-1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(2598.076211353316*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiUx[1]+(750.0*phiUy[0]-750.0*phiUx[0]-1000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2Sq[0]*phiUy[1]-866.0254037844386*rdxCp2Sq[0]*phiUx[1]+((-750.0*phiUy[0])+750.0*phiUx[0]-2600.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[3])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(727.4613391789284*rdxCp2Sq[1]+4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4382.08854314926*rdxCp2[0]*rdxCp2[1])-727.4613391789284*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-201600.0*rdxCp2R3[1])-1259760.0*rdxCp2[0]*rdxCp2Sq[1]-252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(2182.384017536785*rdxCp2[0]*rdxCp2Sq[1]+16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-43647.6803507357*rdxCp2R3[1])-269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]-36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(88200.0*phiC[0]-37800.0*phiLy[0])*rdxCp2R3[1]+(36373.06695894642*rdxCp2[0]*phiUx[1]-16211.99555884469*rdxCp2[0]*phiLy[1]+((-31500.0*phiUx[0])-233370.0*phiLy[0]+642810.0*phiC[0]-252000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(269472.4646415659*rdxCp2Sq[0]*phiUx[1]-2182.384017536785*rdxCp2Sq[0]*phiLy[1]+((-233370.0*phiUx[0])-31500.0*phiLy[0]+642810.0*phiC[0]-1259760.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+43647.6803507357*rdxCp2R3[0]*phiUx[1]+((-37800.0*phiUx[0])+88200.0*phiC[0]-201600.0*bcVals[0])*rdxCp2R3[0])*omega-88200.0*phiC[0]*rdxCp2R3[1]-642810.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-642810.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-88200.0*phiC[0]*rdxCp2R3[0]))/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]-7274.613391789284*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-403221.4280020346*rdxCp2[0]*rdxCp2Sq[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(12600.0*rdxCp2[0]*rdxCp2Sq[1]+8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-93600.0*rdxCp2[0]*rdxCp2Sq[1])-12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(176400.0*phiC[1]-75600.0*phiLy[1])*rdxCp2R3[1]+(210000.0*rdxCp2[0]*phiUx[1]-82980.0*rdxCp2[0]*phiLy[1]+1285620.0*rdxCp2[0]*phiC[1]+((-181865.3347947321*phiUx[0])-81059.97779422344*phiLy[0]+1454922.678357857*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(341400.0*rdxCp2Sq[0]*phiUx[1]-11340.0*rdxCp2Sq[0]*phiLy[1]+1285620.0*rdxCp2Sq[0]*phiC[1]+((-295661.0728520073*phiUx[0])-10911.92008768392*phiLy[0]+1313587.332460236*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+50400.0*rdxCp2R3[0]*phiUx[1]+176400.0*rdxCp2R3[0]*phiC[1]+(174590.7214029428*bcVals[0]-43647.6803507357*phiUx[0])*rdxCp2R3[0])*omega-176400.0*phiC[1]*rdxCp2R3[1]-1285620.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-1285620.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-176400.0*rdxCp2R3[0]*phiC[1]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = (((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+(3780.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-7274.613391789284*rho[0]*rdxCp2Sq[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiUx[3]+((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(174590.7214029428*rdxCp2R3[1]+1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]+1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(11340.0*rdxCp2[0]*rdxCp2Sq[1]+82980.0*rdxCp2Sq[0]*rdxCp2[1]+75600.0*rdxCp2R3[0])*phiUx[2]+((-50400.0*rdxCp2R3[1])-341400.0*rdxCp2[0]*rdxCp2Sq[1]-210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-176400.0*rdxCp2R3[1])-1285620.0*rdxCp2[0]*rdxCp2Sq[1]-1285620.0*rdxCp2Sq[0]*rdxCp2[1]-176400.0*rdxCp2R3[0])*phiC[2]-43647.6803507357*phiLy[0]*rdxCp2R3[1]+(12600.0*rdxCp2[0]*phiUx[1]-8640.0*rdxCp2[0]*phiLy[1]+((-10911.92008768392*phiUx[0])-295661.0728520073*phiLy[0]-87295.36070147139*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(93600.0*rdxCp2Sq[0]*phiUx[1]-12600.0*rdxCp2Sq[0]*phiLy[1]+((-81059.97779422344*phiUx[0])-181865.3347947321*phiLy[0]-403221.4280020346*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0])*phiC[2])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+(640.8587988004846*rdxCp2[0]*rdxCp2[1]+2424.871130596428*rdxCp2Sq[0])*rho[2]+((-2424.871130596428*rdxCp2Sq[1])-640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]-5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiUx[3]+((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-117600.0*rdxCp2R3[1])-857080.0*rdxCp2[0]*rdxCp2Sq[1]-857080.0*rdxCp2Sq[0]*rdxCp2[1]-117600.0*rdxCp2R3[0])*phiC[3]+(168000.0*rdxCp2Sq[0]*rdxCp2[1]-26400.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[3]+(21823.84017536785*rdxCp2[0]*rdxCp2Sq[1]+128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]+29098.45356715714*rdxCp2R3[0])*phiUx[2]+((-16627.68775266122*rdxCp2[0]*rdxCp2Sq[1])-24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]-29098.45356715714*phiLy[1]*rdxCp2R3[1]+(24248.71130596428*rdxCp2[0]*phiUx[1]-128933.8621154272*rdxCp2[0]*phiLy[1]+((-21000.0*phiUx[0])-14400.0*phiLy[0]+168000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(16627.68775266122*rdxCp2Sq[0]*phiUx[1]-21823.84017536785*rdxCp2Sq[0]*phiLy[1]+((-14400.0*phiUx[0])-21000.0*phiLy[0]-26400.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0])*phiC[3])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(207.8460969082653*rdxCp2R3[1]+6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]+13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1143.153532995459*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+1091.192008768392*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(3200.0*rdxCp2R4[1]+96720.0*rdxCp2[0]*rdxCp2R3[1]+222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(623.5382907247956*rdxCp2[0]*rdxCp2R3[1]+22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]+48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(1385.640646055102*rdxCp2R4[1]+42816.29596310264*rdxCp2[0]*rdxCp2R3[1]+122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]+72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(1200.0*phiLy[0]-1200.0*phiC[0])*rdxCp2R4[1]+((-2078.460969082652*rdxCp2[0]*phiUx[1])+2078.460969082652*rdxCp2[0]*phiLy[1]+(1800.0*phiUx[0]+37080.0*phiLy[0]-42480.0*phiC[0]+14400.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-74720.67183852136*rdxCp2Sq[0]*phiUx[1])+6131.459858793824*rdxCp2Sq[0]*phiLy[1]+(64710.0*phiUx[0]+106140.0*phiLy[0]-260670.0*phiC[0]+359280.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-188308.5637988883*rdxCp2R3[0]*phiUx[1])+4364.768035073569*rdxCp2R3[0]*phiLy[1]+(163080.0*phiUx[0]+63000.0*phiLy[0]-446040.0*phiC[0]+879840.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-65471.52052610354*rdxCp2R4[0]*phiUx[1]+(56700.0*phiUx[0]-132300.0*phiC[0]+302400.0*bcVals[0])*rdxCp2R4[0])*omega+1200.0*phiC[0]*rdxCp2R4[1]+42480.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+132300.0*phiC[0]*rdxCp2R4[0])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = (((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(2160.0*rdxCp2[0]*rdxCp2Sq[1]+5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]+11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+10911.92008768392*rdxCp2R3[0]*rho[0])*omega*volFac+(((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+(2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(31869.73485926734*rdxCp2[0]*rdxCp2R3[1]+81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]+14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(3600.0*rdxCp2[0]*rdxCp2R3[1]+25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25920.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(12000.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+25200.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(2400.0*phiLy[1]-2400.0*phiC[1])*rdxCp2R4[1]+((-12000.0*rdxCp2[0]*phiUx[1])+24960.0*rdxCp2[0]*phiLy[1]-84960.0*rdxCp2[0]*phiC[1]+(10392.30484541326*phiUx[0]+10392.30484541326*phiLy[0]-83138.4387633061*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-114600.0*rdxCp2Sq[0]*phiUx[1])+67140.0*rdxCp2Sq[0]*phiLy[1]-521340.0*rdxCp2Sq[0]*phiC[1]+(99246.51127369665*phiUx[0]+30657.29929396912*phiLy[0]-519615.2422706631*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-237600.0*rdxCp2R3[0]*phiUx[1])+22680.0*rdxCp2R3[0]*phiLy[1]-892080.0*rdxCp2R3[0]*phiC[1]+(205767.6359391825*phiUx[0]+21823.84017536785*phiLy[0]-910365.9044582016*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-75600.0*rdxCp2R4[0]*phiUx[1]-264600.0*rdxCp2R4[0]*phiC[1]+(65471.52052610354*phiUx[0]-261886.0821044141*bcVals[0])*rdxCp2R4[0])*omega+2400.0*phiC[1]*rdxCp2R4[1]+84960.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+264600.0*rdxCp2R4[0]*phiC[1])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = (((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+(120.0*rdxCp2R3[1]+3870.0*rdxCp2[0]*rdxCp2Sq[1]+15150.0*rdxCp2Sq[0]*rdxCp2[1]+15750.0*rdxCp2R3[0])*rho[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiUx[3]+(692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(1385.640646055102*rdxCp2R4[1]+43647.6803507357*rdxCp2[0]*rdxCp2R3[1]+145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]+121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(360.0*rdxCp2[0]*rdxCp2R3[1]+12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]+50760.0*rdxCp2R3[0]*rdxCp2[1]+56700.0*rdxCp2R4[0])*phiUx[2]+((-1800.0*rdxCp2[0]*rdxCp2R3[1])-50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-105000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-1200.0*rdxCp2R4[1])-42480.0*rdxCp2[0]*rdxCp2R3[1]-260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]-446040.0*rdxCp2R3[0]*rdxCp2[1]-132300.0*rdxCp2R4[0])*phiC[2]+((-600.0*rdxCp2[0]*phiUx[1])+600.0*rdxCp2[0]*phiLy[1]+(519.6152422706631*phiUx[0]-1558.845726811989*phiLy[0]+4156.921938165305*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((18706.14872174387*phiUx[0]-43647.6803507357*phiLy[0]+99766.1265159673*bcVals[0])*rdxCp2Sq[0]-21600.0*rdxCp2Sq[0]*phiUx[1])*rdxCp2Sq[1]+((-46800.0*rdxCp2R3[0]*phiUx[1])-6300.0*rdxCp2R3[0]*phiLy[1]+(40529.98889711172*phiUx[0]-90932.66739736605*phiLy[0]+201610.7140010173*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0])*phiC[2])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+(727.4613391789284*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiUx[3]+((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2R4[1])-84960.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892080.0*rdxCp2R3[0]*rdxCp2[1]-264600.0*rdxCp2R4[0])*phiC[3]+(10400.0*rdxCp2[0]*rdxCp2R3[1]+35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(2078.460969082652*rdxCp2[0]*rdxCp2R3[1]+39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]+145907.9600296021*rdxCp2R3[0]*rdxCp2[1]+65471.52052610354*rdxCp2R4[0])*phiUx[2]+(3464.101615137754*rdxCp2[0]*rdxCp2R3[1]-36373.06695894642*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-3464.101615137754*rdxCp2[0]*phiUx[1])-17320.50807568877*rdxCp2[0]*phiLy[1]+(3000.0*phiUx[0]+3000.0*phiLy[0]-24000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-24941.53162899183*rdxCp2Sq[0]*phiUx[1])-87295.36070147139*rdxCp2Sq[0]*phiLy[1]+(21600.0*phiUx[0]-86400.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-24941.53162899183*rdxCp2R3[0]*phiUx[1])-32735.76026305177*rdxCp2R3[0]*phiLy[1]+(21600.0*phiUx[0]-31500.0*phiLy[0]+39600.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0])*phiC[3])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-1091.192008768392*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-13146.26562944778*rdxCp2[0]*rdxCp2Sq[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]-207.8460969082653*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+((5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(302400.0*rdxCp2R4[1]+879840.0*rdxCp2[0]*rdxCp2R3[1]+359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-4364.768035073569*rdxCp2[0]*rdxCp2R3[1])-6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(65471.52052610354*rdxCp2R4[1]+188308.5637988883*rdxCp2[0]*rdxCp2R3[1]+74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(56700.0*phiLy[0]-132300.0*phiC[0])*rdxCp2R4[1]+((-72746.13391789283*rdxCp2[0]*phiUx[1])-48635.98667653406*rdxCp2[0]*phiLy[1]+(63000.0*phiUx[0]+163080.0*phiLy[0]-446040.0*phiC[0]-42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-122559.9151435737*rdxCp2Sq[0]*phiUx[1])-22447.37846609264*rdxCp2Sq[0]*phiLy[1]+(106140.0*phiUx[0]+64710.0*phiLy[0]-260670.0*phiC[0]-222560.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-42816.29596310264*rdxCp2R3[0]*phiUx[1])-623.5382907247956*rdxCp2R3[0]*phiLy[1]+(37080.0*phiUx[0]+1800.0*phiLy[0]-42480.0*phiC[0]-96720.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]-1385.640646055102*rdxCp2R4[0]*phiUx[1]+(1200.0*phiUx[0]-1200.0*phiC[0]-3200.0*bcVals[0])*rdxCp2R4[0])*omega+132300.0*phiC[0]*rdxCp2R4[1]+446040.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+1200.0*phiC[0]*rdxCp2R4[0])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-885.0*rdxCp2[0]*rdxCp2Sq[1])-360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-15750.0*rdxCp2R3[1])-15150.0*rdxCp2[0]*rdxCp2Sq[1]-3870.0*rdxCp2Sq[0]*rdxCp2[1]-120.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+((692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiUx[3]+((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]+4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(6300.0*rdxCp2[0]*rdxCp2R3[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(46800.0*rdxCp2[0]*rdxCp2R3[1]+21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(132300.0*phiC[1]-56700.0*phiLy[1])*rdxCp2R4[1]+(105000.0*rdxCp2[0]*phiUx[1]-50760.0*rdxCp2[0]*phiLy[1]+446040.0*rdxCp2[0]*phiC[1]+((-90932.66739736605*phiUx[0])+40529.98889711172*phiLy[0]-121243.5565298214*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+(50400.0*rdxCp2Sq[0]*phiUx[1]-12870.0*rdxCp2Sq[0]*phiLy[1]+260670.0*rdxCp2Sq[0]*phiC[1]+((-43647.6803507357*phiUx[0])+18706.14872174387*phiLy[0]-145838.6779972995*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(1800.0*rdxCp2R3[0]*phiUx[1]-360.0*rdxCp2R3[0]*phiLy[1]+42480.0*rdxCp2R3[0]*phiC[1]+((-1558.845726811989*phiUx[0])+519.6152422706631*phiLy[0]-43647.6803507357*bcVals[0])*rdxCp2R3[0])*rdxCp2[1]+1200.0*rdxCp2R4[0]*phiC[1]-1385.640646055102*bcVals[0]*rdxCp2R4[0])*omega-132300.0*phiC[1]*rdxCp2R4[1]-446040.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-260670.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-42480.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-1200.0*rdxCp2R4[0]*phiC[1]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+((-5670.0*rdxCp2R3[1])-18450.0*rdxCp2[0]*rdxCp2Sq[1]-9480.0*rdxCp2Sq[0]*rdxCp2[1]-1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+10911.92008768392*rho[0]*rdxCp2R3[1]+30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiUx[3]+((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-261886.0821044141*rdxCp2R4[1])-910365.9044582016*rdxCp2[0]*rdxCp2R3[1]-519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]-83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-22680.0*rdxCp2[0]*rdxCp2R3[1])-67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiUx[2]+(75600.0*rdxCp2R4[1]+237600.0*rdxCp2[0]*rdxCp2R3[1]+114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]+12000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[2]+65471.52052610354*phiLy[0]*rdxCp2R4[1]+((-25200.0*rdxCp2[0]*phiUx[1])-25920.0*rdxCp2[0]*phiLy[1]+(21823.84017536785*phiUx[0]+205767.6359391825*phiLy[0]-14549.22678357857*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((-35400.0*rdxCp2Sq[0]*phiUx[1])-25920.0*rdxCp2Sq[0]*phiLy[1]+(30657.29929396912*phiUx[0]+99246.51127369665*phiLy[0]-81752.798117251*bcVals[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-12000.0*rdxCp2R3[0]*phiUx[1])-3600.0*rdxCp2R3[0]*phiLy[1]+(10392.30484541326*phiUx[0]+10392.30484541326*phiLy[0]-31869.73485926734*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[2]))/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+((-961.2881982007268*rdxCp2[0]*rdxCp2Sq[1])-1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiUx[3]+((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[3]+(39600.0*rdxCp2[0]*rdxCp2R3[1]-86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(32735.76026305177*rdxCp2[0]*rdxCp2R3[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]+17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiUx[2]+(24941.53162899183*rdxCp2[0]*rdxCp2R3[1]+24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]-65471.52052610354*phiLy[1]*rdxCp2R4[1]+(36373.06695894642*rdxCp2[0]*phiUx[1]-145907.9600296021*rdxCp2[0]*phiLy[1]+((-31500.0*phiUx[0])+21600.0*phiLy[0]-42000.0*bcVals[0])*rdxCp2[0])*rdxCp2R3[1]+((21600.0*phiLy[0]-35400.0*bcVals[0])*rdxCp2Sq[0]-39386.83536411626*rdxCp2Sq[0]*phiLy[1])*rdxCp2Sq[1]+((-3464.101615137754*rdxCp2R3[0]*phiUx[1])-2078.460969082652*rdxCp2R3[0]*phiLy[1]+(3000.0*phiUx[0]+3000.0*phiLy[0]-10400.0*bcVals[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[3])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-25.98076211353316*rdxCp2Sq[1])-285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+(285.7883832488647*rdxCp2[0]*rdxCp2[1]+25.98076211353316*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+(600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-400.0*rdxCp2R3[1])-4440.0*rdxCp2[0]*rdxCp2Sq[1]-200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-103.9230484541326*rdxCp2[0]*rdxCp2Sq[1])-519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-173.2050807568877*rdxCp2R3[1])-1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(150.0*phiC[0]-150.0*phiLy[0])*rdxCp2R3[1]+(346.4101615137754*rdxCp2[0]*phiUx[1]+519.6152422706631*rdxCp2[0]*phiLy[1]+((-300.0*phiUx[0])-1710.0*phiLy[0]+2010.0*phiC[0]+200.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(1974.53792062852*rdxCp2Sq[0]*phiUx[1]+103.9230484541326*rdxCp2Sq[0]*phiLy[1]+((-1710.0*phiUx[0])-300.0*phiLy[0]+2010.0*phiC[0]+4440.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]+173.2050807568877*rdxCp2R3[0]*phiUx[1]+((-150.0*phiUx[0])+150.0*phiC[0]+400.0*bcVals[0])*rdxCp2R3[0])*omega-150.0*phiC[0]*rdxCp2R3[1]-2010.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-2010.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-150.0*phiC[0]*rdxCp2R3[0]))/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = (((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]-540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]-259.8076211353315*rdxCp2Sq[0]*rho[0])*omega*volFac+((1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiUx[3]+(1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-7967.433714816835*rdxCp2[0]*rdxCp2Sq[1])-346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(900.0*rdxCp2[0]*rdxCp2Sq[1]-900.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(900.0*phiLy[1]-900.0*phiC[1])*rdxCp2R3[1]+((-3000.0*rdxCp2[0]*phiUx[1])+3060.0*rdxCp2[0]*phiLy[1]-12060.0*rdxCp2[0]*phiC[1]+(2598.076211353316*phiUx[0]-2598.076211353316*phiLy[0]+3464.101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+((-600.0*rdxCp2Sq[0]*phiUx[1])+360.0*rdxCp2Sq[0]*phiLy[1]-12060.0*rdxCp2Sq[0]*phiC[1]+(519.6152422706631*phiUx[0]-519.6152422706631*phiLy[0]+12124.35565298214*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1]-900.0*rdxCp2R3[0]*phiC[1]+1039.230484541326*bcVals[0]*rdxCp2R3[0])*omega+900.0*phiC[1]*rdxCp2R3[1]+12060.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+900.0*rdxCp2R3[0]*phiC[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+((-90.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-259.8076211353315*rho[0]*rdxCp2Sq[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiUx[3]+(1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-1039.230484541326*rdxCp2R3[1])-12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-3060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiUx[2]+(600.0*rdxCp2[0]*rdxCp2Sq[1]+3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[2]+(600.0*rdxCp2[0]*phiUx[1]+900.0*rdxCp2[0]*phiLy[1]+((-519.6152422706631*phiUx[0])+519.6152422706631*phiLy[0]+346.4101615137754*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(3000.0*rdxCp2Sq[0]*phiUx[1]-900.0*rdxCp2Sq[0]*phiLy[1]+((-2598.076211353316*phiUx[0])+2598.076211353316*phiLy[0]+7967.433714816835*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[2]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+((-181.8653347947321*rdxCp2[0]*rdxCp2[1])-129.9038105676658*rdxCp2Sq[0])*rho[2]+(129.9038105676658*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUx[3]+((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[3]+((-2600.0*rdxCp2[0]*rdxCp2Sq[1])-1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(519.6152422706631*rdxCp2[0]*rdxCp2Sq[1]+2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiUx[2]+(866.0254037844386*rdxCp2Sq[0]*rdxCp2[1]-866.0254037844386*rdxCp2[0]*rdxCp2Sq[1])*phiLy[2]+((-866.0254037844386*rdxCp2[0]*phiUx[1])-2598.076211353316*rdxCp2[0]*phiLy[1]+(750.0*phiUx[0]-750.0*phiLy[0]+1000.0*bcVals[0])*rdxCp2[0])*rdxCp2Sq[1]+(866.0254037844386*rdxCp2Sq[0]*phiUx[1]-519.6152422706631*rdxCp2Sq[0]*phiLy[1]+((-750.0*phiUx[0])+750.0*phiLy[0]+2600.0*bcVals[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[3])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-727.4613391789284*rdxCp2Sq[1])-4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+(4382.08854314926*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[1]-21000.0*rho[0]*rdxCp2Sq[1]-130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(43647.6803507357*rdxCp2R3[1]+269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-2182.384017536785*rdxCp2[0]*rdxCp2Sq[1])-16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-201600.0*rdxCp2R3[1])-1259760.0*rdxCp2[0]*rdxCp2Sq[1]-252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(88200.0*phiC[0]-37800.0*phiUy[0])*rdxCp2R3[1]+(16211.99555884469*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiLx[1]-252000.0*rdxCp2[0]*bcVals[1]+((-233370.0*phiUy[0])-31500.0*phiLx[0]+642810.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(2182.384017536785*rdxCp2Sq[0]*phiUy[1]-269472.4646415659*rdxCp2Sq[0]*phiLx[1]-1259760.0*rdxCp2Sq[0]*bcVals[1]+((-31500.0*phiUy[0])-233370.0*phiLx[0]+642810.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-43647.6803507357*rdxCp2R3[0]*phiLx[1]-201600.0*rdxCp2R3[0]*bcVals[1]+(88200.0*phiC[0]-37800.0*phiLx[0])*rdxCp2R3[0])*omega-88200.0*phiC[0]*rdxCp2R3[1]-642810.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-642810.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-88200.0*phiC[0]*rdxCp2R3[0]))/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = (((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(21000.0*rdxCp2Sq[1]+26550.0*rdxCp2[0]*rdxCp2[1]+3780.0*rdxCp2Sq[0])*rho[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]-7274.613391789284*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(93600.0*rdxCp2[0]*rdxCp2Sq[1]+12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-12600.0*rdxCp2[0]*rdxCp2Sq[1])-8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-403221.4280020346*rdxCp2[0]*rdxCp2Sq[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(75600.0*phiUy[1]-176400.0*phiC[1])*rdxCp2R3[1]+(82980.0*rdxCp2[0]*phiUy[1]-210000.0*rdxCp2[0]*phiLx[1]-1285620.0*rdxCp2[0]*phiC[1]+1454922.678357857*rdxCp2[0]*bcVals[1]+((-81059.97779422344*phiUy[0])-181865.3347947321*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(11340.0*rdxCp2Sq[0]*phiUy[1]-341400.0*rdxCp2Sq[0]*phiLx[1]-1285620.0*rdxCp2Sq[0]*phiC[1]+1313587.332460236*rdxCp2Sq[0]*bcVals[1]+((-10911.92008768392*phiUy[0])-295661.0728520073*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-50400.0*rdxCp2R3[0]*phiLx[1]-176400.0*rdxCp2R3[0]*phiC[1]+174590.7214029428*rdxCp2R3[0]*bcVals[1]-43647.6803507357*phiLx[0]*rdxCp2R3[0])*omega+176400.0*phiC[1]*rdxCp2R3[1]+1285620.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+176400.0*rdxCp2R3[0]*phiC[1])/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]+1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]-7274.613391789284*rho[0]*rdxCp2Sq[1]-43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiLx[3]+(50400.0*rdxCp2R3[1]+341400.0*rdxCp2[0]*rdxCp2Sq[1]+210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-11340.0*rdxCp2[0]*rdxCp2Sq[1])-82980.0*rdxCp2Sq[0]*rdxCp2[1]-75600.0*rdxCp2R3[0])*phiLx[2]+(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0])*phiC[2]+(174590.7214029428*rdxCp2R3[1]+1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]+1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]-43647.6803507357*phiUy[0]*rdxCp2R3[1]+(8640.0*rdxCp2[0]*phiUy[1]-12600.0*rdxCp2[0]*phiLx[1]-87295.36070147139*rdxCp2[0]*bcVals[1]+((-295661.0728520073*phiUy[0])-10911.92008768392*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(12600.0*rdxCp2Sq[0]*phiUy[1]-93600.0*rdxCp2Sq[0]*phiLx[1]-403221.4280020346*rdxCp2Sq[0]*bcVals[1]+((-181865.3347947321*phiUy[0])-81059.97779422344*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-176400.0*rdxCp2R3[1])-1285620.0*rdxCp2[0]*rdxCp2Sq[1]-1285620.0*rdxCp2Sq[0]*rdxCp2[1]-176400.0*rdxCp2R3[0])*phiC[2]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-640.8587988004846*rdxCp2[0]*rdxCp2[1])-2424.871130596428*rdxCp2Sq[0])*rho[2]+(2424.871130596428*rdxCp2Sq[1]+640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]-5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiLx[3]+((-117600.0*rdxCp2R3[1])-857080.0*rdxCp2[0]*rdxCp2Sq[1]-857080.0*rdxCp2Sq[0]*rdxCp2[1]-117600.0*rdxCp2R3[0])*phiC[3]+(16627.68775266122*rdxCp2[0]*rdxCp2Sq[1]+24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-21823.84017536785*rdxCp2[0]*rdxCp2Sq[1])-128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]-29098.45356715714*rdxCp2R3[0])*phiLx[2]+(168000.0*rdxCp2Sq[0]*rdxCp2[1]-26400.0*rdxCp2[0]*rdxCp2Sq[1])*bcVals[2]+29098.45356715714*phiUy[1]*rdxCp2R3[1]+(128933.8621154272*rdxCp2[0]*phiUy[1]-24248.71130596428*rdxCp2[0]*phiLx[1]+168000.0*rdxCp2[0]*bcVals[1]+((-14400.0*phiUy[0])-21000.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(21823.84017536785*rdxCp2Sq[0]*phiUy[1]-16627.68775266122*rdxCp2Sq[0]*phiLx[1]-26400.0*rdxCp2Sq[0]*bcVals[1]+((-21000.0*phiUy[0])-14400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0])*phiC[3])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-207.8460969082653*rdxCp2R3[1])-6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]-13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1143.153532995459*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-1091.192008768392*rdxCp2R3[0])*rho[1]+1200.0*rho[0]*rdxCp2R3[1]+36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+((2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-1385.640646055102*rdxCp2R4[1])-42816.29596310264*rdxCp2[0]*rdxCp2R3[1]-122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]-72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-623.5382907247956*rdxCp2[0]*rdxCp2R3[1])-22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]-48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-3200.0*rdxCp2R4[1])-96720.0*rdxCp2[0]*rdxCp2R3[1]-222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(1200.0*phiUy[0]-1200.0*phiC[0])*rdxCp2R4[1]+((-2078.460969082652*rdxCp2[0]*phiUy[1])+2078.460969082652*rdxCp2[0]*phiLx[1]+14400.0*rdxCp2[0]*bcVals[1]+(37080.0*phiUy[0]+1800.0*phiLx[0]-42480.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-6131.459858793824*rdxCp2Sq[0]*phiUy[1])+74720.67183852136*rdxCp2Sq[0]*phiLx[1]+359280.0*rdxCp2Sq[0]*bcVals[1]+(106140.0*phiUy[0]+64710.0*phiLx[0]-260670.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-4364.768035073569*rdxCp2R3[0]*phiUy[1])+188308.5637988883*rdxCp2R3[0]*phiLx[1]+879840.0*rdxCp2R3[0]*bcVals[1]+(63000.0*phiUy[0]+163080.0*phiLx[0]-446040.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+65471.52052610354*rdxCp2R4[0]*phiLx[1]+302400.0*rdxCp2R4[0]*bcVals[1]+(56700.0*phiLx[0]-132300.0*phiC[0])*rdxCp2R4[0])*omega+1200.0*phiC[0]*rdxCp2R4[1]+42480.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+132300.0*phiC[0]*rdxCp2R4[0])/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-1200.0*rdxCp2R3[1])-9480.0*rdxCp2[0]*rdxCp2Sq[1]-18450.0*rdxCp2Sq[0]*rdxCp2[1]-5670.0*rdxCp2R3[0])*rho[1]+11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+10911.92008768392*rdxCp2R3[0]*rho[0])*omega*volFac+((2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-12000.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25200.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-3600.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+((-31869.73485926734*rdxCp2[0]*rdxCp2R3[1])-81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]-14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(2400.0*phiC[1]-2400.0*phiUy[1])*rdxCp2R4[1]+((-24960.0*rdxCp2[0]*phiUy[1])+12000.0*rdxCp2[0]*phiLx[1]+84960.0*rdxCp2[0]*phiC[1]-83138.4387633061*rdxCp2[0]*bcVals[1]+(10392.30484541326*phiUy[0]+10392.30484541326*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-67140.0*rdxCp2Sq[0]*phiUy[1])+114600.0*rdxCp2Sq[0]*phiLx[1]+521340.0*rdxCp2Sq[0]*phiC[1]-519615.2422706631*rdxCp2Sq[0]*bcVals[1]+(30657.29929396912*phiUy[0]+99246.51127369665*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-22680.0*rdxCp2R3[0]*phiUy[1])+237600.0*rdxCp2R3[0]*phiLx[1]+892080.0*rdxCp2R3[0]*phiC[1]-910365.9044582016*rdxCp2R3[0]*bcVals[1]+(21823.84017536785*phiUy[0]+205767.6359391825*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]+75600.0*rdxCp2R4[0]*phiLx[1]+264600.0*rdxCp2R4[0]*phiC[1]-261886.0821044141*rdxCp2R4[0]*bcVals[1]+65471.52052610354*phiLx[0]*rdxCp2R4[0])*omega-2400.0*phiC[1]*rdxCp2R4[1]-84960.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-892080.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-264600.0*rdxCp2R4[0]*phiC[1]))/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+((-120.0*rdxCp2R3[1])-3870.0*rdxCp2[0]*rdxCp2Sq[1]-15150.0*rdxCp2Sq[0]*rdxCp2[1]-15750.0*rdxCp2R3[0])*rho[2]+((-360.0*rdxCp2[0]*rdxCp2Sq[1])-885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+346.4101615137754*rho[0]*rdxCp2R3[1]+10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[3]+(1800.0*rdxCp2[0]*rdxCp2R3[1]+50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+105000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]-50760.0*rdxCp2R3[0]*rdxCp2[1]-56700.0*rdxCp2R4[0])*phiLx[2]+(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0])*phiC[2]+((-1385.640646055102*rdxCp2R4[1])-43647.6803507357*rdxCp2[0]*rdxCp2R3[1]-145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]-121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+((-600.0*rdxCp2[0]*phiUy[1])+600.0*rdxCp2[0]*phiLx[1]+4156.921938165305*rdxCp2[0]*bcVals[1]+(519.6152422706631*phiLx[0]-1558.845726811989*phiUy[0])*rdxCp2[0])*rdxCp2R3[1]+(21600.0*rdxCp2Sq[0]*phiLx[1]+99766.1265159673*rdxCp2Sq[0]*bcVals[1]+(18706.14872174387*phiLx[0]-43647.6803507357*phiUy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(6300.0*rdxCp2R3[0]*phiUy[1]+46800.0*rdxCp2R3[0]*phiLx[1]+201610.7140010173*rdxCp2R3[0]*bcVals[1]+(40529.98889711172*phiLx[0]-90932.66739736605*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-1200.0*rdxCp2R4[1])-42480.0*rdxCp2[0]*rdxCp2R3[1]-260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]-446040.0*rdxCp2R3[0]*rdxCp2[1]-132300.0*rdxCp2R4[0])*phiC[2]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+((-346.4101615137754*rdxCp2R3[1])-1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]-961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiLx[3]+((-2400.0*rdxCp2R4[1])-84960.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892080.0*rdxCp2R3[0]*rdxCp2[1]-264600.0*rdxCp2R4[0])*phiC[3]+(36373.06695894642*rdxCp2R3[0]*rdxCp2[1]-3464.101615137754*rdxCp2[0]*rdxCp2R3[1])*phiUy[2]+((-2078.460969082652*rdxCp2[0]*rdxCp2R3[1])-39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]-145907.9600296021*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[2]+((-10400.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(17320.50807568877*rdxCp2[0]*phiUy[1]+3464.101615137754*rdxCp2[0]*phiLx[1]-24000.0*rdxCp2[0]*bcVals[1]+(3000.0*phiUy[0]+3000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(87295.36070147139*rdxCp2Sq[0]*phiUy[1]+24941.53162899183*rdxCp2Sq[0]*phiLx[1]-86400.0*rdxCp2Sq[0]*bcVals[1]+21600.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(32735.76026305177*rdxCp2R3[0]*phiUy[1]+24941.53162899183*rdxCp2R3[0]*phiLx[1]+39600.0*rdxCp2R3[0]*bcVals[1]+(21600.0*phiLx[0]-31500.0*phiUy[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0])*phiC[3])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = (((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1091.192008768392*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(13146.26562944778*rdxCp2[0]*rdxCp2Sq[1]+6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[1]+31500.0*rho[0]*rdxCp2R3[1]+91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-65471.52052610354*rdxCp2R4[1])-188308.5637988883*rdxCp2[0]*rdxCp2R3[1]-74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(4364.768035073569*rdxCp2[0]*rdxCp2R3[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(302400.0*rdxCp2R4[1]+879840.0*rdxCp2[0]*rdxCp2R3[1]+359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]+14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(56700.0*phiUy[0]-132300.0*phiC[0])*rdxCp2R4[1]+(48635.98667653406*rdxCp2[0]*phiUy[1]+72746.13391789283*rdxCp2[0]*phiLx[1]+42000.0*rdxCp2[0]*bcVals[1]+(163080.0*phiUy[0]+63000.0*phiLx[0]-446040.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(22447.37846609264*rdxCp2Sq[0]*phiUy[1]+122559.9151435737*rdxCp2Sq[0]*phiLx[1]+222560.0*rdxCp2Sq[0]*bcVals[1]+(64710.0*phiUy[0]+106140.0*phiLx[0]-260670.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(623.5382907247956*rdxCp2R3[0]*phiUy[1]+42816.29596310264*rdxCp2R3[0]*phiLx[1]+96720.0*rdxCp2R3[0]*bcVals[1]+(1800.0*phiUy[0]+37080.0*phiLx[0]-42480.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]+1385.640646055102*rdxCp2R4[0]*phiLx[1]+3200.0*rdxCp2R4[0]*bcVals[1]+(1200.0*phiLx[0]-1200.0*phiC[0])*rdxCp2R4[0])*omega+132300.0*phiC[0]*rdxCp2R4[1]+446040.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]+260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]+1200.0*phiC[0]*rdxCp2R4[0])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = (((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(885.0*rdxCp2[0]*rdxCp2Sq[1]+360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(15750.0*rdxCp2R3[1]+15150.0*rdxCp2[0]*rdxCp2Sq[1]+3870.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[1]+21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+(((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-46800.0*rdxCp2[0]*rdxCp2R3[1])-21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(600.0*rdxCp2R3[0]*rdxCp2[1]-6300.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(201610.7140010173*rdxCp2[0]*rdxCp2R3[1]+99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]+4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+(56700.0*phiUy[1]-132300.0*phiC[1])*rdxCp2R4[1]+(50760.0*rdxCp2[0]*phiUy[1]-105000.0*rdxCp2[0]*phiLx[1]-446040.0*rdxCp2[0]*phiC[1]+121243.5565298214*rdxCp2[0]*bcVals[1]+(40529.98889711172*phiUy[0]-90932.66739736605*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(12870.0*rdxCp2Sq[0]*phiUy[1]-50400.0*rdxCp2Sq[0]*phiLx[1]-260670.0*rdxCp2Sq[0]*phiC[1]+145838.6779972995*rdxCp2Sq[0]*bcVals[1]+(18706.14872174387*phiUy[0]-43647.6803507357*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(360.0*rdxCp2R3[0]*phiUy[1]-1800.0*rdxCp2R3[0]*phiLx[1]-42480.0*rdxCp2R3[0]*phiC[1]+43647.6803507357*rdxCp2R3[0]*bcVals[1]+(519.6152422706631*phiUy[0]-1558.845726811989*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-1200.0*rdxCp2R4[0]*phiC[1]+1385.640646055102*rdxCp2R4[0]*bcVals[1])*omega+132300.0*phiC[1]*rdxCp2R4[1]+446040.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+1200.0*rdxCp2R4[0]*phiC[1])/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = (((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+(5310.0*rdxCp2[0]*rdxCp2Sq[1]+2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+10911.92008768392*rho[0]*rdxCp2R3[1]+30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiLx[3]+((-75600.0*rdxCp2R4[1])-237600.0*rdxCp2[0]*rdxCp2R3[1]-114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-12000.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+(22680.0*rdxCp2[0]*rdxCp2R3[1]+67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiLx[2]+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[2]+((-261886.0821044141*rdxCp2R4[1])-910365.9044582016*rdxCp2[0]*rdxCp2R3[1]-519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]-83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+65471.52052610354*phiUy[0]*rdxCp2R4[1]+(25920.0*rdxCp2[0]*phiUy[1]+25200.0*rdxCp2[0]*phiLx[1]+14549.22678357857*rdxCp2[0]*bcVals[1]+(205767.6359391825*phiUy[0]+21823.84017536785*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(25920.0*rdxCp2Sq[0]*phiUy[1]+35400.0*rdxCp2Sq[0]*phiLx[1]+81752.798117251*rdxCp2Sq[0]*bcVals[1]+(99246.51127369665*phiUy[0]+30657.29929396912*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(3600.0*rdxCp2R3[0]*phiUy[1]+12000.0*rdxCp2R3[0]*phiLx[1]+31869.73485926734*rdxCp2R3[0]*bcVals[1]+(10392.30484541326*phiUy[0]+10392.30484541326*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[2])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+(5455.960043841962*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]+8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]+3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiUy[3]+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[3]+((-24941.53162899183*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiUy[2]+((-32735.76026305177*rdxCp2[0]*rdxCp2R3[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]-17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(39600.0*rdxCp2[0]*rdxCp2R3[1]-86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[2]+65471.52052610354*phiUy[1]*rdxCp2R4[1]+(145907.9600296021*rdxCp2[0]*phiUy[1]-36373.06695894642*rdxCp2[0]*phiLx[1]+42000.0*rdxCp2[0]*bcVals[1]+(21600.0*phiUy[0]-31500.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(39386.83536411626*rdxCp2Sq[0]*phiUy[1]+35400.0*rdxCp2Sq[0]*bcVals[1]+21600.0*phiUy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+(2078.460969082652*rdxCp2R3[0]*phiUy[1]+3464.101615137754*rdxCp2R3[0]*phiLx[1]+10400.0*rdxCp2R3[0]*bcVals[1]+(3000.0*phiUy[0]+3000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[3])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(25.98076211353316*rdxCp2Sq[1]+285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+((-285.7883832488647*rdxCp2[0]*rdxCp2[1])-25.98076211353316*rdxCp2Sq[0])*rho[1]-150.0*rho[0]*rdxCp2Sq[1]-1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]-150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(173.2050807568877*rdxCp2R3[1]+1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(400.0*rdxCp2R3[1]+4440.0*rdxCp2[0]*rdxCp2Sq[1]+200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(150.0*phiC[0]-150.0*phiUy[0])*rdxCp2R3[1]+((-519.6152422706631*rdxCp2[0]*phiUy[1])-346.4101615137754*rdxCp2[0]*phiLx[1]-200.0*rdxCp2[0]*bcVals[1]+((-1710.0*phiUy[0])-300.0*phiLx[0]+2010.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-103.9230484541326*rdxCp2Sq[0]*phiUy[1])-1974.53792062852*rdxCp2Sq[0]*phiLx[1]-4440.0*rdxCp2Sq[0]*bcVals[1]+((-300.0*phiUy[0])-1710.0*phiLx[0]+2010.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]-173.2050807568877*rdxCp2R3[0]*phiLx[1]-400.0*rdxCp2R3[0]*bcVals[1]+(150.0*phiC[0]-150.0*phiLx[0])*rdxCp2R3[0])*omega-150.0*phiC[0]*rdxCp2R3[1]-2010.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]-2010.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]-150.0*phiC[0]*rdxCp2R3[0]))/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-450.0*rdxCp2Sq[1])-1080.0*rdxCp2[0]*rdxCp2[1]-90.0*rdxCp2Sq[0])*rho[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]-259.8076211353315*rdxCp2Sq[0]*rho[0])*omega*volFac+((1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(3000.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(900.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(7967.433714816835*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(900.0*phiC[1]-900.0*phiUy[1])*rdxCp2R3[1]+((-3060.0*rdxCp2[0]*phiUy[1])+3000.0*rdxCp2[0]*phiLx[1]+12060.0*rdxCp2[0]*phiC[1]-3464.101615137754*rdxCp2[0]*bcVals[1]+(2598.076211353316*phiLx[0]-2598.076211353316*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-360.0*rdxCp2Sq[0]*phiUy[1])+600.0*rdxCp2Sq[0]*phiLx[1]+12060.0*rdxCp2Sq[0]*phiC[1]-12124.35565298214*rdxCp2Sq[0]*bcVals[1]+(519.6152422706631*phiLx[0]-519.6152422706631*phiUy[0])*rdxCp2Sq[0])*rdxCp2[1]+900.0*rdxCp2R3[0]*phiC[1]-1039.230484541326*rdxCp2R3[0]*bcVals[1])*omega-900.0*phiC[1]*rdxCp2R3[1]-12060.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-900.0*rdxCp2R3[0]*phiC[1]))/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = (((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]-540.0*rdxCp2[0]*rdxCp2[1]*rho[1]-259.8076211353315*rho[0]*rdxCp2Sq[1]-2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiLx[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+3060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiLx[2]+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[2]+(1039.230484541326*rdxCp2R3[1]+12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+((-900.0*rdxCp2[0]*phiUy[1])-600.0*rdxCp2[0]*phiLx[1]-346.4101615137754*rdxCp2[0]*bcVals[1]+(519.6152422706631*phiUy[0]-519.6152422706631*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(900.0*rdxCp2Sq[0]*phiUy[1]-3000.0*rdxCp2Sq[0]*phiLx[1]-7967.433714816835*rdxCp2Sq[0]*bcVals[1]+(2598.076211353316*phiUy[0]-2598.076211353316*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[2])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0])*rho[2]+((-129.9038105676658*rdxCp2Sq[1])-181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]-900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiUy[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[3]+(866.0254037844386*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*phiUy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2Sq[1])-2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(2600.0*rdxCp2[0]*rdxCp2Sq[1]+1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[2]+(2598.076211353316*rdxCp2[0]*phiUy[1]+866.0254037844386*rdxCp2[0]*phiLx[1]-1000.0*rdxCp2[0]*bcVals[1]+(750.0*phiLx[0]-750.0*phiUy[0])*rdxCp2[0])*rdxCp2Sq[1]+(519.6152422706631*rdxCp2Sq[0]*phiUy[1]-866.0254037844386*rdxCp2Sq[0]*phiLx[1]-2600.0*rdxCp2Sq[0]*bcVals[1]+(750.0*phiUy[0]-750.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[3])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((177.0*rdxCp2[0]*rdxCp2[1]*rho[3]+((-727.4613391789284*rdxCp2Sq[1])-4382.08854314926*rdxCp2[0]*rdxCp2[1])*rho[2]+((-4382.08854314926*rdxCp2[0]*rdxCp2[1])-727.4613391789284*rdxCp2Sq[0])*rho[1]+21000.0*rho[0]*rdxCp2Sq[1]+130280.0*rdxCp2[0]*rho[0]*rdxCp2[1]+21000.0*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-18720.0*rdxCp2[0]*rdxCp2Sq[1])-2520.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-2520.0*rdxCp2[0]*rdxCp2Sq[1])-18720.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(201600.0*rdxCp2R3[1]+1259760.0*rdxCp2[0]*rdxCp2Sq[1]+252000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(43647.6803507357*rdxCp2R3[1]+269472.4646415659*rdxCp2[0]*rdxCp2Sq[1]+36373.06695894642*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-2182.384017536785*rdxCp2[0]*rdxCp2Sq[1])-16211.99555884469*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(37800.0*phiLy[0]-88200.0*phiC[0])*rdxCp2R3[1]+((-16211.99555884469*rdxCp2[0]*phiLy[1])+36373.06695894642*rdxCp2[0]*phiLx[1]+252000.0*rdxCp2[0]*bcVals[1]+(233370.0*phiLy[0]+31500.0*phiLx[0]-642810.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+((-2182.384017536785*rdxCp2Sq[0]*phiLy[1])+269472.4646415659*rdxCp2Sq[0]*phiLx[1]+1259760.0*rdxCp2Sq[0]*bcVals[1]+(31500.0*phiLy[0]+233370.0*phiLx[0]-642810.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+43647.6803507357*rdxCp2R3[0]*phiLx[1]+201600.0*rdxCp2R3[0]*bcVals[1]+(37800.0*phiLx[0]-88200.0*phiC[0])*rdxCp2R3[0])*omega+88200.0*phiC[0]*rdxCp2R3[1]+642810.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+642810.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+88200.0*phiC[0]*rdxCp2R3[0])/(88200.0*rdxCp2R3[1]+642810.0*rdxCp2[0]*rdxCp2Sq[1]+642810.0*rdxCp2Sq[0]*rdxCp2[1]+88200.0*rdxCp2R3[0]); 
  phiC[1] = -(1.0*(((727.4613391789284*rdxCp2Sq[1]+192.2576396401454*rdxCp2[0]*rdxCp2[1])*rho[3]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[2]+((-21000.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-3780.0*rdxCp2Sq[0])*rho[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1]+7274.613391789284*rdxCp2Sq[0]*rho[0])*omega*volFac+(((-87295.36070147139*rdxCp2R3[1])-95817.05067471028*rdxCp2[0]*rdxCp2Sq[1]-13094.30410522071*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-14549.22678357857*rdxCp2[0]*rdxCp2Sq[1])-9976.61265159673*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(403221.4280020346*rdxCp2[0]*rdxCp2Sq[1]+87295.36070147139*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(93600.0*rdxCp2[0]*rdxCp2Sq[1]+12600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-12600.0*rdxCp2[0]*rdxCp2Sq[1])-8640.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(176400.0*phiC[1]-75600.0*phiLy[1])*rdxCp2R3[1]+((-82980.0*rdxCp2[0]*phiLy[1])+210000.0*rdxCp2[0]*phiLx[1]+1285620.0*rdxCp2[0]*phiC[1]-1454922.678357857*rdxCp2[0]*bcVals[1]+(81059.97779422344*phiLy[0]+181865.3347947321*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-11340.0*rdxCp2Sq[0]*phiLy[1])+341400.0*rdxCp2Sq[0]*phiLx[1]+1285620.0*rdxCp2Sq[0]*phiC[1]-1313587.332460236*rdxCp2Sq[0]*bcVals[1]+(10911.92008768392*phiLy[0]+295661.0728520073*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]+50400.0*rdxCp2R3[0]*phiLx[1]+176400.0*rdxCp2R3[0]*phiC[1]-174590.7214029428*rdxCp2R3[0]*bcVals[1]+43647.6803507357*phiLx[0]*rdxCp2R3[0])*omega-176400.0*phiC[1]*rdxCp2R3[1]-1285620.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]-1285620.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]-176400.0*rdxCp2R3[0]*phiC[1]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[2] = -(1.0*(((192.2576396401454*rdxCp2[0]*rdxCp2[1]+727.4613391789284*rdxCp2Sq[0])*rho[3]+((-3780.0*rdxCp2Sq[1])-26550.0*rdxCp2[0]*rdxCp2[1]-21000.0*rdxCp2Sq[0])*rho[2]-1770.0*rdxCp2[0]*rdxCp2[1]*rho[1]+7274.613391789284*rho[0]*rdxCp2Sq[1]+43820.88543149259*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-9976.61265159673*rdxCp2[0]*rdxCp2Sq[1])-14549.22678357857*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-13094.30410522071*rdxCp2[0]*rdxCp2Sq[1])-95817.05067471028*rdxCp2Sq[0]*rdxCp2[1]-87295.36070147139*rdxCp2R3[0])*phiLx[3]+((-174590.7214029428*rdxCp2R3[1])-1313587.332460236*rdxCp2[0]*rdxCp2Sq[1]-1454922.678357857*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(50400.0*rdxCp2R3[1]+341400.0*rdxCp2[0]*rdxCp2Sq[1]+210000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-11340.0*rdxCp2[0]*rdxCp2Sq[1])-82980.0*rdxCp2Sq[0]*rdxCp2[1]-75600.0*rdxCp2R3[0])*phiLx[2]+(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0])*phiC[2]+43647.6803507357*phiLy[0]*rdxCp2R3[1]+((-8640.0*rdxCp2[0]*phiLy[1])+12600.0*rdxCp2[0]*phiLx[1]+87295.36070147139*rdxCp2[0]*bcVals[1]+(295661.0728520073*phiLy[0]+10911.92008768392*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-12600.0*rdxCp2Sq[0]*phiLy[1])+93600.0*rdxCp2Sq[0]*phiLx[1]+403221.4280020346*rdxCp2Sq[0]*bcVals[1]+(181865.3347947321*phiLy[0]+81059.97779422344*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+((-176400.0*rdxCp2R3[1])-1285620.0*rdxCp2[0]*rdxCp2Sq[1]-1285620.0*rdxCp2Sq[0]*rdxCp2[1]-176400.0*rdxCp2R3[0])*phiC[2]))/(176400.0*rdxCp2R3[1]+1285620.0*rdxCp2[0]*rdxCp2Sq[1]+1285620.0*rdxCp2Sq[0]*rdxCp2[1]+176400.0*rdxCp2R3[0]); 
  phiC[3] = (((1260.0*rdxCp2Sq[1]+7333.0*rdxCp2[0]*rdxCp2[1]+1260.0*rdxCp2Sq[0])*rho[3]+((-640.8587988004846*rdxCp2[0]*rdxCp2[1])-2424.871130596428*rdxCp2Sq[0])*rho[2]+((-2424.871130596428*rdxCp2Sq[1])-640.8587988004846*rdxCp2[0]*rdxCp2[1])*rho[1]+5900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-33600.0*rdxCp2R3[1])-148880.0*rdxCp2[0]*rdxCp2Sq[1]-25200.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-25200.0*rdxCp2[0]*rdxCp2Sq[1])-148880.0*rdxCp2Sq[0]*rdxCp2[1]-33600.0*rdxCp2R3[0])*phiLx[3]+((-117600.0*rdxCp2R3[1])-857080.0*rdxCp2[0]*rdxCp2Sq[1]-857080.0*rdxCp2Sq[0]*rdxCp2[1]-117600.0*rdxCp2R3[0])*phiC[3]+(26400.0*rdxCp2[0]*rdxCp2Sq[1]-168000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(16627.68775266122*rdxCp2[0]*rdxCp2Sq[1]+24248.71130596428*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-21823.84017536785*rdxCp2[0]*rdxCp2Sq[1])-128933.8621154272*rdxCp2Sq[0]*rdxCp2[1]-29098.45356715714*rdxCp2R3[0])*phiLx[2]-29098.45356715714*phiLy[1]*rdxCp2R3[1]+((-128933.8621154272*rdxCp2[0]*phiLy[1])+24248.71130596428*rdxCp2[0]*phiLx[1]-168000.0*rdxCp2[0]*bcVals[1]+(14400.0*phiLy[0]+21000.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-21823.84017536785*rdxCp2Sq[0]*phiLy[1])+16627.68775266122*rdxCp2Sq[0]*phiLx[1]+26400.0*rdxCp2Sq[0]*bcVals[1]+(21000.0*phiLy[0]+14400.0*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0])*phiC[3])/(117600.0*rdxCp2R3[1]+857080.0*rdxCp2[0]*rdxCp2Sq[1]+857080.0*rdxCp2Sq[0]*rdxCp2[1]+117600.0*rdxCp2R3[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((216.0*rdxCp2[0]*rdxCp2Sq[1]+531.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-207.8460969082653*rdxCp2R3[1])-6235.382907247957*rdxCp2[0]*rdxCp2Sq[1]-13146.26562944778*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1143.153532995459*rdxCp2[0]*rdxCp2Sq[1]+3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]+1091.192008768392*rdxCp2R3[0])*rho[1]-1200.0*rho[0]*rdxCp2R3[1]-36540.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-91020.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-31500.0*rdxCp2R3[0]*rho[0])*omega*volFac+((2400.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+5040.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-720.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-56160.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-3200.0*rdxCp2R4[1])-96720.0*rdxCp2[0]*rdxCp2R3[1]-222560.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-1385.640646055102*rdxCp2R4[1])-42816.29596310264*rdxCp2[0]*rdxCp2R3[1]-122559.9151435737*rdxCp2Sq[0]*rdxCp2Sq[1]-72746.13391789283*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-623.5382907247956*rdxCp2[0]*rdxCp2R3[1])-22447.37846609264*rdxCp2Sq[0]*rdxCp2Sq[1]-48635.98667653406*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(1200.0*phiC[0]-1200.0*phiLy[0])*rdxCp2R4[1]+(2078.460969082652*rdxCp2[0]*phiLy[1]-2078.460969082652*rdxCp2[0]*phiLx[1]-14400.0*rdxCp2[0]*bcVals[1]+((-37080.0*phiLy[0])-1800.0*phiLx[0]+42480.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+(6131.459858793824*rdxCp2Sq[0]*phiLy[1]-74720.67183852136*rdxCp2Sq[0]*phiLx[1]-359280.0*rdxCp2Sq[0]*bcVals[1]+((-106140.0*phiLy[0])-64710.0*phiLx[0]+260670.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(4364.768035073569*rdxCp2R3[0]*phiLy[1]-188308.5637988883*rdxCp2R3[0]*phiLx[1]-879840.0*rdxCp2R3[0]*bcVals[1]+((-63000.0*phiLy[0])-163080.0*phiLx[0]+446040.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-65471.52052610354*rdxCp2R4[0]*phiLx[1]-302400.0*rdxCp2R4[0]*bcVals[1]+(132300.0*phiC[0]-56700.0*phiLx[0])*rdxCp2R4[0])*omega-1200.0*phiC[0]*rdxCp2R4[1]-42480.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-446040.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-132300.0*phiC[0]*rdxCp2R4[0]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[1] = (((207.8460969082653*rdxCp2R3[1]+1122.368923304632*rdxCp2[0]*rdxCp2Sq[1]+576.7729189204359*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+((-2160.0*rdxCp2[0]*rdxCp2Sq[1])-5310.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+(1200.0*rdxCp2R3[1]+9480.0*rdxCp2[0]*rdxCp2Sq[1]+18450.0*rdxCp2Sq[0]*rdxCp2[1]+5670.0*rdxCp2R3[0])*rho[1]-11431.53532995459*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-30657.29929396912*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-10911.92008768392*rdxCp2R3[0]*rho[0])*omega*volFac+((2771.281292110204*rdxCp2R4[1]+28821.32543794612*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+26188.60821044141*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-4156.921938165305*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-29929.83795479019*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-31869.73485926734*rdxCp2[0]*rdxCp2R3[1])-81752.798117251*rdxCp2Sq[0]*rdxCp2Sq[1]-14549.22678357857*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-12000.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25200.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-3600.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-25920.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(2400.0*phiLy[1]-2400.0*phiC[1])*rdxCp2R4[1]+(24960.0*rdxCp2[0]*phiLy[1]-12000.0*rdxCp2[0]*phiLx[1]-84960.0*rdxCp2[0]*phiC[1]+83138.4387633061*rdxCp2[0]*bcVals[1]+((-10392.30484541326*phiLy[0])-10392.30484541326*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+(67140.0*rdxCp2Sq[0]*phiLy[1]-114600.0*rdxCp2Sq[0]*phiLx[1]-521340.0*rdxCp2Sq[0]*phiC[1]+519615.2422706631*rdxCp2Sq[0]*bcVals[1]+((-30657.29929396912*phiLy[0])-99246.51127369665*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+(22680.0*rdxCp2R3[0]*phiLy[1]-237600.0*rdxCp2R3[0]*phiLx[1]-892080.0*rdxCp2R3[0]*phiC[1]+910365.9044582016*rdxCp2R3[0]*bcVals[1]+((-21823.84017536785*phiLy[0])-205767.6359391825*phiLx[0])*rdxCp2R3[0])*rdxCp2[1]-75600.0*rdxCp2R4[0]*phiLx[1]-264600.0*rdxCp2R4[0]*phiC[1]+261886.0821044141*rdxCp2R4[0]*bcVals[1]-65471.52052610354*phiLx[0]*rdxCp2R4[0])*omega+2400.0*phiC[1]*rdxCp2R4[1]+84960.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]+264600.0*rdxCp2R4[0]*phiC[1])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 
  phiC[2] = -(1.0*(((72.74613391789283*rdxCp2[0]*rdxCp2Sq[1]+306.5729929396912*rdxCp2Sq[0]*rdxCp2[1]+545.5960043841961*rdxCp2R3[0])*rho[3]+((-120.0*rdxCp2R3[1])-3870.0*rdxCp2[0]*rdxCp2Sq[1]-15150.0*rdxCp2Sq[0]*rdxCp2[1]-15750.0*rdxCp2R3[0])*rho[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+885.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-346.4101615137754*rho[0]*rdxCp2R3[1]-10392.30484541326*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-21910.4427157463*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+((692.8203230275509*rdxCp2[0]*rdxCp2R3[1]-7274.613391789284*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-415.6921938165305*rdxCp2[0]*rdxCp2R3[1])-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-58612.59932813079*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[3]+((-1385.640646055102*rdxCp2R4[1])-43647.6803507357*rdxCp2[0]*rdxCp2R3[1]-145838.6779972995*rdxCp2Sq[0]*rdxCp2Sq[1]-121243.5565298214*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(1800.0*rdxCp2[0]*rdxCp2R3[1]+50400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+105000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-360.0*rdxCp2[0]*rdxCp2R3[1])-12870.0*rdxCp2Sq[0]*rdxCp2Sq[1]-50760.0*rdxCp2R3[0]*rdxCp2[1]-56700.0*rdxCp2R4[0])*phiLx[2]+(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0])*phiC[2]+(600.0*rdxCp2[0]*phiLy[1]-600.0*rdxCp2[0]*phiLx[1]-4156.921938165305*rdxCp2[0]*bcVals[1]+(1558.845726811989*phiLy[0]-519.6152422706631*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-21600.0*rdxCp2Sq[0]*phiLx[1])-99766.1265159673*rdxCp2Sq[0]*bcVals[1]+(43647.6803507357*phiLy[0]-18706.14872174387*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-6300.0*rdxCp2R3[0]*phiLy[1])-46800.0*rdxCp2R3[0]*phiLx[1]-201610.7140010173*rdxCp2R3[0]*bcVals[1]+(90932.66739736605*phiLy[0]-40529.98889711172*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+((-1200.0*rdxCp2R4[1])-42480.0*rdxCp2[0]*rdxCp2R3[1]-260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]-446040.0*rdxCp2R3[0]*rdxCp2[1]-132300.0*rdxCp2R4[0])*phiC[2]))/(1200.0*rdxCp2R4[1]+42480.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+446040.0*rdxCp2R3[0]*rdxCp2[1]+132300.0*rdxCp2R4[0]); 
  phiC[3] = (((120.0*rdxCp2R3[1]+2148.0*rdxCp2[0]*rdxCp2Sq[1]+7893.0*rdxCp2Sq[0]*rdxCp2[1]+2835.0*rdxCp2R3[0])*rho[3]+((-727.4613391789284*rdxCp2[0]*rdxCp2Sq[1])-3065.729929396912*rdxCp2Sq[0]*rdxCp2[1]-5455.960043841962*rdxCp2R3[0])*rho[2]+(346.4101615137754*rdxCp2R3[1]+1870.614872174387*rdxCp2[0]*rdxCp2Sq[1]+961.2881982007268*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-3600.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-8850.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-20000.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-37800.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-2400.0*rdxCp2[0]*rdxCp2R3[1])-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-168480.0*rdxCp2R3[0]*rdxCp2[1]-75600.0*rdxCp2R4[0])*phiLx[3]+((-2400.0*rdxCp2R4[1])-84960.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-892080.0*rdxCp2R3[0]*rdxCp2[1]-264600.0*rdxCp2R4[0])*phiC[3]+((-10400.0*rdxCp2[0]*rdxCp2R3[1])-35400.0*rdxCp2Sq[0]*rdxCp2Sq[1]-42000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+(36373.06695894642*rdxCp2R3[0]*rdxCp2[1]-3464.101615137754*rdxCp2[0]*rdxCp2R3[1])*phiLy[2]+((-2078.460969082652*rdxCp2[0]*rdxCp2R3[1])-39386.83536411626*rdxCp2Sq[0]*rdxCp2Sq[1]-145907.9600296021*rdxCp2R3[0]*rdxCp2[1]-65471.52052610354*rdxCp2R4[0])*phiLx[2]+((-17320.50807568877*rdxCp2[0]*phiLy[1])-3464.101615137754*rdxCp2[0]*phiLx[1]+24000.0*rdxCp2[0]*bcVals[1]+((-3000.0*phiLy[0])-3000.0*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-87295.36070147139*rdxCp2Sq[0]*phiLy[1])-24941.53162899183*rdxCp2Sq[0]*phiLx[1]+86400.0*rdxCp2Sq[0]*bcVals[1]-21600.0*phiLx[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-32735.76026305177*rdxCp2R3[0]*phiLy[1])-24941.53162899183*rdxCp2R3[0]*phiLx[1]-39600.0*rdxCp2R3[0]*bcVals[1]+(31500.0*phiLy[0]-21600.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0])*phiC[3])/(2400.0*rdxCp2R4[1]+84960.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+892080.0*rdxCp2R3[0]*rdxCp2[1]+264600.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = -(1.0*(((531.0*rdxCp2[0]*rdxCp2Sq[1]+216.0*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(1091.192008768392*rdxCp2R3[1]+3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]+1143.153532995459*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-13146.26562944778*rdxCp2[0]*rdxCp2Sq[1])-6235.382907247957*rdxCp2Sq[0]*rdxCp2[1]-207.8460969082653*rdxCp2R3[0])*rho[1]-31500.0*rho[0]*rdxCp2R3[1]-91020.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-36540.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-1200.0*rdxCp2R3[0]*rho[0])*omega*volFac+(((-56160.0*rdxCp2[0]*rdxCp2R3[1])-25920.0*rdxCp2Sq[0]*rdxCp2Sq[1]-720.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(5040.0*rdxCp2[0]*rdxCp2R3[1]+7080.0*rdxCp2Sq[0]*rdxCp2Sq[1]+2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-302400.0*rdxCp2R4[1])-879840.0*rdxCp2[0]*rdxCp2R3[1]-359280.0*rdxCp2Sq[0]*rdxCp2Sq[1]-14400.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-65471.52052610354*rdxCp2R4[1])-188308.5637988883*rdxCp2[0]*rdxCp2R3[1]-74720.67183852136*rdxCp2Sq[0]*rdxCp2Sq[1]-2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(4364.768035073569*rdxCp2[0]*rdxCp2R3[1]+6131.459858793824*rdxCp2Sq[0]*rdxCp2Sq[1]+2078.460969082652*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]+(132300.0*phiC[0]-56700.0*phiLy[0])*rdxCp2R4[1]+((-48635.98667653406*rdxCp2[0]*phiLy[1])-72746.13391789283*rdxCp2[0]*phiLx[1]-42000.0*rdxCp2[0]*bcVals[1]+((-163080.0*phiLy[0])-63000.0*phiLx[0]+446040.0*phiC[0])*rdxCp2[0])*rdxCp2R3[1]+((-22447.37846609264*rdxCp2Sq[0]*phiLy[1])-122559.9151435737*rdxCp2Sq[0]*phiLx[1]-222560.0*rdxCp2Sq[0]*bcVals[1]+((-64710.0*phiLy[0])-106140.0*phiLx[0]+260670.0*phiC[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-623.5382907247956*rdxCp2R3[0]*phiLy[1])-42816.29596310264*rdxCp2R3[0]*phiLx[1]-96720.0*rdxCp2R3[0]*bcVals[1]+((-1800.0*phiLy[0])-37080.0*phiLx[0]+42480.0*phiC[0])*rdxCp2R3[0])*rdxCp2[1]-1385.640646055102*rdxCp2R4[0]*phiLx[1]-3200.0*rdxCp2R4[0]*bcVals[1]+(1200.0*phiC[0]-1200.0*phiLx[0])*rdxCp2R4[0])*omega-132300.0*phiC[0]*rdxCp2R4[1]-446040.0*phiC[0]*rdxCp2[0]*rdxCp2R3[1]-260670.0*phiC[0]*rdxCp2Sq[0]*rdxCp2Sq[1]-42480.0*phiC[0]*rdxCp2R3[0]*rdxCp2[1]-1200.0*phiC[0]*rdxCp2R4[0]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[1] = -(1.0*(((545.5960043841961*rdxCp2R3[1]+306.5729929396912*rdxCp2[0]*rdxCp2Sq[1]+72.74613391789283*rdxCp2Sq[0]*rdxCp2[1])*rho[3]+(885.0*rdxCp2[0]*rdxCp2Sq[1]+360.0*rdxCp2Sq[0]*rdxCp2[1])*rho[2]+((-15750.0*rdxCp2R3[1])-15150.0*rdxCp2[0]*rdxCp2Sq[1]-3870.0*rdxCp2Sq[0]*rdxCp2[1]-120.0*rdxCp2R3[0])*rho[1]-21910.4427157463*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-10392.30484541326*rdxCp2Sq[0]*rho[0]*rdxCp2[1]-346.4101615137754*rdxCp2R3[0]*rho[0])*omega*volFac+(((-65471.52052610354*rdxCp2R4[1])-58612.59932813079*rdxCp2[0]*rdxCp2R3[1]-14860.99592894097*rdxCp2Sq[0]*rdxCp2Sq[1]-415.6921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(692.8203230275509*rdxCp2R3[0]*rdxCp2[1]-7274.613391789284*rdxCp2[0]*rdxCp2R3[1])*phiLx[3]+((-201610.7140010173*rdxCp2[0]*rdxCp2R3[1])-99766.1265159673*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-46800.0*rdxCp2[0]*rdxCp2R3[1])-21600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-600.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(600.0*rdxCp2R3[0]*rdxCp2[1]-6300.0*rdxCp2[0]*rdxCp2R3[1])*phiLx[2]+(132300.0*phiC[1]-56700.0*phiLy[1])*rdxCp2R4[1]+((-50760.0*rdxCp2[0]*phiLy[1])+105000.0*rdxCp2[0]*phiLx[1]+446040.0*rdxCp2[0]*phiC[1]-121243.5565298214*rdxCp2[0]*bcVals[1]+(90932.66739736605*phiLx[0]-40529.98889711172*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-12870.0*rdxCp2Sq[0]*phiLy[1])+50400.0*rdxCp2Sq[0]*phiLx[1]+260670.0*rdxCp2Sq[0]*phiC[1]-145838.6779972995*rdxCp2Sq[0]*bcVals[1]+(43647.6803507357*phiLx[0]-18706.14872174387*phiLy[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-360.0*rdxCp2R3[0]*phiLy[1])+1800.0*rdxCp2R3[0]*phiLx[1]+42480.0*rdxCp2R3[0]*phiC[1]-43647.6803507357*rdxCp2R3[0]*bcVals[1]+(1558.845726811989*phiLx[0]-519.6152422706631*phiLy[0])*rdxCp2R3[0])*rdxCp2[1]+1200.0*rdxCp2R4[0]*phiC[1]-1385.640646055102*rdxCp2R4[0]*bcVals[1])*omega-132300.0*phiC[1]*rdxCp2R4[1]-446040.0*rdxCp2[0]*phiC[1]*rdxCp2R3[1]-260670.0*rdxCp2Sq[0]*phiC[1]*rdxCp2Sq[1]-42480.0*rdxCp2R3[0]*phiC[1]*rdxCp2[1]-1200.0*rdxCp2R4[0]*phiC[1]))/(132300.0*rdxCp2R4[1]+446040.0*rdxCp2[0]*rdxCp2R3[1]+260670.0*rdxCp2Sq[0]*rdxCp2Sq[1]+42480.0*rdxCp2R3[0]*rdxCp2[1]+1200.0*rdxCp2R4[0]); 
  phiC[2] = (((576.7729189204359*rdxCp2[0]*rdxCp2Sq[1]+1122.368923304632*rdxCp2Sq[0]*rdxCp2[1]+207.8460969082653*rdxCp2R3[0])*rho[3]+(5670.0*rdxCp2R3[1]+18450.0*rdxCp2[0]*rdxCp2Sq[1]+9480.0*rdxCp2Sq[0]*rdxCp2[1]+1200.0*rdxCp2R3[0])*rho[2]+((-5310.0*rdxCp2[0]*rdxCp2Sq[1])-2160.0*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-10911.92008768392*rho[0]*rdxCp2R3[1]-30657.29929396912*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-11431.53532995459*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-29929.83795479019*rdxCp2[0]*rdxCp2R3[1])-29929.83795479019*rdxCp2Sq[0]*rdxCp2Sq[1]-4156.921938165305*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+(26188.60821044141*rdxCp2[0]*rdxCp2R3[1]+77526.59414678294*rdxCp2Sq[0]*rdxCp2Sq[1]+28821.32543794612*rdxCp2R3[0]*rdxCp2[1]+2771.281292110204*rdxCp2R4[0])*phiLx[3]+(261886.0821044141*rdxCp2R4[1]+910365.9044582016*rdxCp2[0]*rdxCp2R3[1]+519615.2422706631*rdxCp2Sq[0]*rdxCp2Sq[1]+83138.4387633061*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-75600.0*rdxCp2R4[1])-237600.0*rdxCp2[0]*rdxCp2R3[1]-114600.0*rdxCp2Sq[0]*rdxCp2Sq[1]-12000.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+(22680.0*rdxCp2[0]*rdxCp2R3[1]+67140.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiLx[2]+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[2]-65471.52052610354*phiLy[0]*rdxCp2R4[1]+((-25920.0*rdxCp2[0]*phiLy[1])-25200.0*rdxCp2[0]*phiLx[1]-14549.22678357857*rdxCp2[0]*bcVals[1]+((-205767.6359391825*phiLy[0])-21823.84017536785*phiLx[0])*rdxCp2[0])*rdxCp2R3[1]+((-25920.0*rdxCp2Sq[0]*phiLy[1])-35400.0*rdxCp2Sq[0]*phiLx[1]-81752.798117251*rdxCp2Sq[0]*bcVals[1]+((-99246.51127369665*phiLy[0])-30657.29929396912*phiLx[0])*rdxCp2Sq[0])*rdxCp2Sq[1]+((-3600.0*rdxCp2R3[0]*phiLy[1])-12000.0*rdxCp2R3[0]*phiLx[1]-31869.73485926734*rdxCp2R3[0]*bcVals[1]+((-10392.30484541326*phiLy[0])-10392.30484541326*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[2])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 
  phiC[3] = (((2835.0*rdxCp2R3[1]+7893.0*rdxCp2[0]*rdxCp2Sq[1]+2148.0*rdxCp2Sq[0]*rdxCp2[1]+120.0*rdxCp2R3[0])*rho[3]+(961.2881982007268*rdxCp2[0]*rdxCp2Sq[1]+1870.614872174387*rdxCp2Sq[0]*rdxCp2[1]+346.4101615137754*rdxCp2R3[0])*rho[2]+((-5455.960043841962*rdxCp2R3[1])-3065.729929396912*rdxCp2[0]*rdxCp2Sq[1]-727.4613391789284*rdxCp2Sq[0]*rdxCp2[1])*rho[1]-8850.0*rdxCp2[0]*rho[0]*rdxCp2Sq[1]-3600.0*rdxCp2Sq[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-75600.0*rdxCp2R4[1])-168480.0*rdxCp2[0]*rdxCp2R3[1]-45480.0*rdxCp2Sq[0]*rdxCp2Sq[1]-2400.0*rdxCp2R3[0]*rdxCp2[1])*phiLy[3]+((-37800.0*rdxCp2[0]*rdxCp2R3[1])-100800.0*rdxCp2Sq[0]*rdxCp2Sq[1]-20000.0*rdxCp2R3[0]*rdxCp2[1])*phiLx[3]+((-264600.0*rdxCp2R4[1])-892080.0*rdxCp2[0]*rdxCp2R3[1]-521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]-84960.0*rdxCp2R3[0]*rdxCp2[1]-2400.0*rdxCp2R4[0])*phiC[3]+((-39600.0*rdxCp2[0]*rdxCp2R3[1])+86400.0*rdxCp2Sq[0]*rdxCp2Sq[1]+24000.0*rdxCp2R3[0]*rdxCp2[1])*bcVals[3]+((-24941.53162899183*rdxCp2[0]*rdxCp2R3[1])-24941.53162899183*rdxCp2Sq[0]*rdxCp2Sq[1]-3464.101615137754*rdxCp2R3[0]*rdxCp2[1])*phiLy[2]+((-32735.76026305177*rdxCp2[0]*rdxCp2R3[1])-87295.36070147139*rdxCp2Sq[0]*rdxCp2Sq[1]-17320.50807568877*rdxCp2R3[0]*rdxCp2[1])*phiLx[2]-65471.52052610354*phiLy[1]*rdxCp2R4[1]+((-145907.9600296021*rdxCp2[0]*phiLy[1])+36373.06695894642*rdxCp2[0]*phiLx[1]-42000.0*rdxCp2[0]*bcVals[1]+(31500.0*phiLx[0]-21600.0*phiLy[0])*rdxCp2[0])*rdxCp2R3[1]+((-39386.83536411626*rdxCp2Sq[0]*phiLy[1])-35400.0*rdxCp2Sq[0]*bcVals[1]-21600.0*phiLy[0]*rdxCp2Sq[0])*rdxCp2Sq[1]+((-2078.460969082652*rdxCp2R3[0]*phiLy[1])-3464.101615137754*rdxCp2R3[0]*phiLx[1]-10400.0*rdxCp2R3[0]*bcVals[1]+((-3000.0*phiLy[0])-3000.0*phiLx[0])*rdxCp2R3[0])*rdxCp2[1])*omega+(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0])*phiC[3])/(264600.0*rdxCp2R4[1]+892080.0*rdxCp2[0]*rdxCp2R3[1]+521340.0*rdxCp2Sq[0]*rdxCp2Sq[1]+84960.0*rdxCp2R3[0]*rdxCp2[1]+2400.0*rdxCp2R4[0]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  phiC[0] = ((54.0*rdxCp2[0]*rdxCp2[1]*rho[3]+(25.98076211353316*rdxCp2Sq[1]+285.7883832488647*rdxCp2[0]*rdxCp2[1])*rho[2]+(285.7883832488647*rdxCp2[0]*rdxCp2[1]+25.98076211353316*rdxCp2Sq[0])*rho[1]+150.0*rho[0]*rdxCp2Sq[1]+1680.0*rdxCp2[0]*rho[0]*rdxCp2[1]+150.0*rdxCp2Sq[0]*rho[0])*omega*volFac+((600.0*rdxCp2[0]*rdxCp2Sq[1]+120.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(120.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+(400.0*rdxCp2R3[1]+4440.0*rdxCp2[0]*rdxCp2Sq[1]+200.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(173.2050807568877*rdxCp2R3[1]+1974.53792062852*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(103.9230484541326*rdxCp2[0]*rdxCp2Sq[1]+519.6152422706631*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+(150.0*phiLy[0]-150.0*phiC[0])*rdxCp2R3[1]+(519.6152422706631*rdxCp2[0]*phiLy[1]+346.4101615137754*rdxCp2[0]*phiLx[1]+200.0*rdxCp2[0]*bcVals[1]+(1710.0*phiLy[0]+300.0*phiLx[0]-2010.0*phiC[0])*rdxCp2[0])*rdxCp2Sq[1]+(103.9230484541326*rdxCp2Sq[0]*phiLy[1]+1974.53792062852*rdxCp2Sq[0]*phiLx[1]+4440.0*rdxCp2Sq[0]*bcVals[1]+(300.0*phiLy[0]+1710.0*phiLx[0]-2010.0*phiC[0])*rdxCp2Sq[0])*rdxCp2[1]+173.2050807568877*rdxCp2R3[0]*phiLx[1]+400.0*rdxCp2R3[0]*bcVals[1]+(150.0*phiLx[0]-150.0*phiC[0])*rdxCp2R3[0])*omega+150.0*phiC[0]*rdxCp2R3[1]+2010.0*phiC[0]*rdxCp2[0]*rdxCp2Sq[1]+2010.0*phiC[0]*rdxCp2Sq[0]*rdxCp2[1]+150.0*phiC[0]*rdxCp2R3[0])/(150.0*rdxCp2R3[1]+2010.0*rdxCp2[0]*rdxCp2Sq[1]+2010.0*rdxCp2Sq[0]*rdxCp2[1]+150.0*rdxCp2R3[0]); 
  phiC[1] = (((77.94228634059945*rdxCp2Sq[1]+109.1192008768392*rdxCp2[0]*rdxCp2[1])*rho[3]+540.0*rdxCp2[0]*rdxCp2[1]*rho[2]+(450.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+90.0*rdxCp2Sq[0])*rho[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1]+259.8076211353315*rdxCp2Sq[0]*rho[0])*omega*volFac+((1039.230484541326*rdxCp2R3[1]+3533.383647440509*rdxCp2[0]*rdxCp2Sq[1]+415.6921938165305*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(1039.230484541326*rdxCp2Sq[0]*rdxCp2[1]-1039.230484541326*rdxCp2[0]*rdxCp2Sq[1])*phiLx[3]+(7967.433714816835*rdxCp2[0]*rdxCp2Sq[1]+346.4101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(3000.0*rdxCp2[0]*rdxCp2Sq[1]+600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(900.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2[0]*rdxCp2Sq[1])*phiLx[2]+(900.0*phiLy[1]-900.0*phiC[1])*rdxCp2R3[1]+(3060.0*rdxCp2[0]*phiLy[1]-3000.0*rdxCp2[0]*phiLx[1]-12060.0*rdxCp2[0]*phiC[1]+3464.101615137754*rdxCp2[0]*bcVals[1]+(2598.076211353316*phiLy[0]-2598.076211353316*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+(360.0*rdxCp2Sq[0]*phiLy[1]-600.0*rdxCp2Sq[0]*phiLx[1]-12060.0*rdxCp2Sq[0]*phiC[1]+12124.35565298214*rdxCp2Sq[0]*bcVals[1]+(519.6152422706631*phiLy[0]-519.6152422706631*phiLx[0])*rdxCp2Sq[0])*rdxCp2[1]-900.0*rdxCp2R3[0]*phiC[1]+1039.230484541326*rdxCp2R3[0]*bcVals[1])*omega+900.0*phiC[1]*rdxCp2R3[1]+12060.0*rdxCp2[0]*phiC[1]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*phiC[1]*rdxCp2[1]+900.0*rdxCp2R3[0]*phiC[1])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[2] = (((109.1192008768392*rdxCp2[0]*rdxCp2[1]+77.94228634059945*rdxCp2Sq[0])*rho[3]+(90.0*rdxCp2Sq[1]+1080.0*rdxCp2[0]*rdxCp2[1]+450.0*rdxCp2Sq[0])*rho[2]+540.0*rdxCp2[0]*rdxCp2[1]*rho[1]+259.8076211353315*rho[0]*rdxCp2Sq[1]+2857.883832488647*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+((1039.230484541326*rdxCp2[0]*rdxCp2Sq[1]-1039.230484541326*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+(415.6921938165305*rdxCp2[0]*rdxCp2Sq[1]+3533.383647440509*rdxCp2Sq[0]*rdxCp2[1]+1039.230484541326*rdxCp2R3[0])*phiLx[3]+(1039.230484541326*rdxCp2R3[1]+12124.35565298214*rdxCp2[0]*rdxCp2Sq[1]+3464.101615137754*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+(360.0*rdxCp2[0]*rdxCp2Sq[1]+3060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiLx[2]+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[2]+(900.0*rdxCp2[0]*phiLy[1]+600.0*rdxCp2[0]*phiLx[1]+346.4101615137754*rdxCp2[0]*bcVals[1]+(519.6152422706631*phiLx[0]-519.6152422706631*phiLy[0])*rdxCp2[0])*rdxCp2Sq[1]+((-900.0*rdxCp2Sq[0]*phiLy[1])+3000.0*rdxCp2Sq[0]*phiLx[1]+7967.433714816835*rdxCp2Sq[0]*bcVals[1]+(2598.076211353316*phiLx[0]-2598.076211353316*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[2])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 
  phiC[3] = (((45.0*rdxCp2Sq[1]+288.0*rdxCp2[0]*rdxCp2[1]+45.0*rdxCp2Sq[0])*rho[3]+(181.8653347947321*rdxCp2[0]*rdxCp2[1]+129.9038105676658*rdxCp2Sq[0])*rho[2]+(129.9038105676658*rdxCp2Sq[1]+181.8653347947321*rdxCp2[0]*rdxCp2[1])*rho[1]+900.0*rdxCp2[0]*rho[0]*rdxCp2[1])*omega*volFac+(((-3000.0*rdxCp2[0]*rdxCp2Sq[1])-600.0*rdxCp2Sq[0]*rdxCp2[1])*phiLy[3]+((-600.0*rdxCp2[0]*rdxCp2Sq[1])-3000.0*rdxCp2Sq[0]*rdxCp2[1])*phiLx[3]+((-900.0*rdxCp2R3[1])-12060.0*rdxCp2[0]*rdxCp2Sq[1]-12060.0*rdxCp2Sq[0]*rdxCp2[1]-900.0*rdxCp2R3[0])*phiC[3]+(2600.0*rdxCp2[0]*rdxCp2Sq[1]+1000.0*rdxCp2Sq[0]*rdxCp2[1])*bcVals[3]+(866.0254037844386*rdxCp2[0]*rdxCp2Sq[1]-866.0254037844386*rdxCp2Sq[0]*rdxCp2[1])*phiLy[2]+((-519.6152422706631*rdxCp2[0]*rdxCp2Sq[1])-2598.076211353316*rdxCp2Sq[0]*rdxCp2[1])*phiLx[2]+((-2598.076211353316*rdxCp2[0]*phiLy[1])-866.0254037844386*rdxCp2[0]*phiLx[1]+1000.0*rdxCp2[0]*bcVals[1]+(750.0*phiLy[0]-750.0*phiLx[0])*rdxCp2[0])*rdxCp2Sq[1]+((-519.6152422706631*rdxCp2Sq[0]*phiLy[1])+866.0254037844386*rdxCp2Sq[0]*phiLx[1]+2600.0*rdxCp2Sq[0]*bcVals[1]+(750.0*phiLx[0]-750.0*phiLy[0])*rdxCp2Sq[0])*rdxCp2[1])*omega+(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0])*phiC[3])/(900.0*rdxCp2R3[1]+12060.0*rdxCp2[0]*rdxCp2Sq[1]+12060.0*rdxCp2Sq[0]*rdxCp2[1]+900.0*rdxCp2R3[0]); 

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
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 
  rdxLxCu[0] = rdxLxSq[0]*rdxLx[0]; 
  rdxUxCu[0] = rdxUxSq[0]*rdxUx[0]; 
  rdxLxR4[0] = rdxLxCu[0]*rdxLx[0]; 
  rdxUxR4[0] = rdxUxCu[0]*rdxUx[0]; 
  rdxLx[1]   = 1.0/dxLx[1]; 
  rdxUx[1]   = 1.0/dxUx[1]; 
  rdxLxSq[1] = rdxLx[1]*rdxLx[1]; 
  rdxUxSq[1] = rdxUx[1]*rdxUx[1]; 
  rdxLxCu[1] = rdxLxSq[1]*rdxLx[1]; 
  rdxUxCu[1] = rdxUxSq[1]*rdxUx[1]; 
  rdxLxR4[1] = rdxLxCu[1]*rdxLx[1]; 
  rdxUxR4[1] = rdxUxCu[1]*rdxUx[1]; 
  rdxCp2[1]  = volFac*4.0/(dxC[1]*dxC[1]); 
  rdxLy[0]   = 1.0/dxLy[0]; 
  rdxUy[0]   = 1.0/dxUy[0]; 
  rdxLySq[0] = rdxLy[0]*rdxLy[0]; 
  rdxUySq[0] = rdxUy[0]*rdxUy[0]; 
  rdxLyCu[0] = rdxLySq[0]*rdxLy[0]; 
  rdxUyCu[0] = rdxUySq[0]*rdxUy[0]; 
  rdxLyR4[0] = rdxLyCu[0]*rdxLy[0]; 
  rdxUyR4[0] = rdxUyCu[0]*rdxUy[0]; 
  rdxLy[1]   = 1.0/dxLy[1]; 
  rdxUy[1]   = 1.0/dxUy[1]; 
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
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-17.32050807568877*rdxUy[1]*phiUy[3]+17.32050807568877*rdxLy[1]*phiLy[3]+(17.32050807568877*rdxLy[1]-17.32050807568877*rdxUy[1])*phiC[3]+(18.0*phiUy[1]-18.0*phiC[1])*rdxUy[1]+(18.0*phiLy[1]-18.0*phiC[1])*rdxLy[1]-7.0*rdxUx[0]*phiUx[1]-7.0*rdxLx[0]*phiLx[1]+((-23.0*rdxUx[0])-23.0*rdxLx[0])*phiC[1]+(8.660254037844386*phiUx[0]-22.5166604983954*phiC[0])*rdxUx[0]+(22.5166604983954*phiC[0]-8.660254037844386*phiLx[0])*rdxLx[0]); 
  resOut[2] = 0.0625*(16.0*rho[2]*volFac-17.32050807568877*rdxUx[0]*phiUx[3]+17.32050807568877*rdxLx[0]*phiLx[3]+(17.32050807568877*rdxLx[0]-17.32050807568877*rdxUx[0])*phiC[3]-7.0*rdxUy[1]*phiUy[2]+18.0*rdxUx[0]*phiUx[2]-7.0*rdxLy[1]*phiLy[2]+18.0*rdxLx[0]*phiLx[2]+((-23.0*rdxUy[1])-23.0*rdxLy[1]-18.0*rdxUx[0]-18.0*rdxLx[0])*phiC[2]+(8.660254037844386*phiUy[0]-22.5166604983954*phiC[0])*rdxUy[1]+(22.5166604983954*phiC[0]-8.660254037844386*phiLy[0])*rdxLy[1]); 
  resOut[3] = 0.125*(8.0*rho[3]*volFac-7.0*rdxUy[1]*phiUy[3]-7.0*rdxUx[0]*phiUx[3]-7.0*rdxLy[1]*phiLy[3]-7.0*rdxLx[0]*phiLx[3]+((-23.0*rdxUy[1])-23.0*rdxLy[1]-23.0*rdxUx[0]-23.0*rdxLx[0])*phiC[3]+8.660254037844386*rdxUx[0]*phiUx[2]-8.660254037844386*rdxLx[0]*phiLx[2]+(22.5166604983954*rdxLx[0]-22.5166604983954*rdxUx[0])*phiC[2]+(8.660254037844386*phiUy[1]-22.5166604983954*phiC[1])*rdxUy[1]+(22.5166604983954*phiC[1]-8.660254037844386*phiLy[1])*rdxLy[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiLy[2]+(3.0*phiUy[0]+3.0*phiLy[0]-6.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-9.0*phiC[0]+24.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiLy[3]+(3.0*phiUy[1]+3.0*phiLy[1]-6.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-50.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]+8.660254037844386*phiC[0]-69.28203230275508*bcVals[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiUx[2]-10.0*rdxCp2[1]*phiLy[2]+((-40.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiUx[3]-20.0*rdxCp2[1]*phiLy[3]+((-80.0*rdxCp2[1])-100.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]+17.32050807568877*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiLy[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]-13.85640646055102*rdxCp2[0]*phiUx[1]-20.78460969082652*rdxCp2[0]*phiC[1]+(12.0*phiUx[0]-12.0*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiLy[3]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.3333333333333333*(3.0*rho[2]*volFac-13.85640646055102*rdxCp2[0]*phiUx[3]-20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+12.0*rdxCp2[0]*phiUx[2]-30.0*rdxCp2[1]*phiLy[2]+((-120.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+(25.98076211353316*phiUy[0]-25.98076211353316*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[1]*phiUy[3]-40.0*rdxCp2[0]*phiUx[3]-60.0*rdxCp2[1]*phiLy[3]+((-240.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+34.64101615137754*rdxCp2[0]*phiUx[2]-34.64101615137754*rdxCp2[0]*phiC[2]+(51.96152422706631*phiUy[1]-51.96152422706631*phiLy[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiLy[2]+(3.0*phiUy[0]+3.0*phiLy[0]-6.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]-3.464101615137754*rdxCp2[0]*phiC[1]+24.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiLy[3]+(3.0*phiUy[1]+3.0*phiLy[1]-6.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiLx[1]-50.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+((-8.660254037844386*phiLx[0])-8.660254037844386*phiC[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac+3.464101615137754*rdxCp2[0]*phiLx[3]-3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]-10.0*rdxCp2[1]*phiLy[2]+3.0*rdxCp2[0]*phiLx[2]+((-40.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+(8.660254037844386*phiUy[0]-8.660254037844386*phiLy[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[1]*phiLy[3]-20.0*rdxCp2[0]*phiLx[3]+((-80.0*rdxCp2[1])-100.0*rdxCp2[0])*phiC[3]-17.32050807568877*rdxCp2[0]*phiLx[2]-17.32050807568877*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]-17.32050807568877*phiLy[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiLy[2]+(9.0*phiUy[0]+9.0*phiLy[0]-18.0*phiC[0])*rdxCp2[1]+13.85640646055102*rdxCp2[0]*phiLx[1]+20.78460969082652*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(12.0*phiLx[0]-12.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiLy[3]+(27.0*phiUy[1]+27.0*phiLy[1]-54.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.3333333333333333*(3.0*rho[2]*volFac+13.85640646055102*rdxCp2[0]*phiLx[3]+20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]-30.0*rdxCp2[1]*phiLy[2]+12.0*rdxCp2[0]*phiLx[2]+((-120.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+(25.98076211353316*phiUy[0]-25.98076211353316*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[1]*phiUy[3]-60.0*rdxCp2[1]*phiLy[3]-40.0*rdxCp2[0]*phiLx[3]+((-240.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]-34.64101615137754*rdxCp2[0]*phiLx[2]+34.64101615137754*rdxCp2[0]*phiC[2]+(51.96152422706631*phiUy[1]-51.96152422706631*phiLy[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiC[2]+24.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiLx[1]+(3.0*phiUx[0]+3.0*phiLx[0]-6.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-10.0*rdxCp2[0]*phiLx[1]-40.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiLx[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiUx[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-6.0*rdxCp2[0])*phiC[2]-69.28203230275508*rdxCp2[1]*bcVals[2]+(8.660254037844386*phiUy[0]+8.660254037844386*phiC[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiUx[3]-20.0*rdxCp2[0]*phiLx[3]+((-100.0*rdxCp2[1])-80.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiLx[2]+(17.32050807568877*phiUy[1]+17.32050807568877*phiC[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-13.85640646055102*rdxCp2[1]*phiUy[2]-20.78460969082652*rdxCp2[1]*phiC[2]-8.0*rdxCp2[1]*bcVals[2]+(12.0*phiUy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.3333333333333333*(3.0*rho[1]*volFac-13.85640646055102*rdxCp2[1]*phiUy[3]-20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiUy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]-120.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiLx[3]-60.0*rdxCp2[1]*phiUy[2]+27.0*rdxCp2[0]*phiUx[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiUy[3]-60.0*rdxCp2[0]*phiUx[3]-60.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-240.0*rdxCp2[0])*phiC[3]+51.96152422706631*rdxCp2[0]*phiUx[2]-51.96152422706631*rdxCp2[0]*phiLx[2]+(34.64101615137754*phiUy[1]-34.64101615137754*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac+24.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]-3.464101615137754*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiLx[1]+(3.0*phiUx[0]+3.0*phiLx[0]-6.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac+3.464101615137754*rdxCp2[1]*phiLy[3]-3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-10.0*rdxCp2[0]*phiLx[1]-40.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiLx[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+3.0*rdxCp2[0]*phiUx[2]-10.0*rdxCp2[1]*phiLy[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-6.0*rdxCp2[0])*phiC[2]+((-8.660254037844386*phiLy[0])-8.660254037844386*phiC[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[0]*phiUx[3]-20.0*rdxCp2[1]*phiLy[3]-20.0*rdxCp2[0]*phiLx[3]+((-100.0*rdxCp2[1])-80.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]-17.32050807568877*rdxCp2[0]*phiLx[2]+((-17.32050807568877*phiLy[1])-17.32050807568877*phiC[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac+8.0*rdxCp2[1]*bcVals[3]+13.85640646055102*rdxCp2[1]*phiLy[2]+20.78460969082652*rdxCp2[1]*phiC[2]+(12.0*phiLy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.3333333333333333*(3.0*rho[1]*volFac+13.85640646055102*rdxCp2[1]*phiLy[3]+20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiLy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-30.0*rdxCp2[0]*phiLx[1]-120.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]-25.98076211353316*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiLx[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+27.0*rdxCp2[0]*phiUx[2]-60.0*rdxCp2[1]*phiLy[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-54.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[0]*phiUx[3]-40.0*rdxCp2[1]*phiLy[3]-60.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-240.0*rdxCp2[0])*phiC[3]+51.96152422706631*rdxCp2[0]*phiUx[2]-51.96152422706631*rdxCp2[0]*phiLx[2]+(34.64101615137754*phiC[1]-34.64101615137754*phiLy[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiC[2]+24.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-9.0*phiC[0]+24.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-50.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]+8.660254037844386*phiC[0]-69.28203230275508*bcVals[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiUx[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]-69.28203230275508*rdxCp2[1]*bcVals[2]+(8.660254037844386*phiUy[0]+8.660254037844386*phiC[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiUx[3]+((-100.0*rdxCp2[1])-100.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]+17.32050807568877*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]+17.32050807568877*phiC[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-13.85640646055102*rdxCp2[1]*phiUy[2]-20.78460969082652*rdxCp2[1]*phiC[2]-8.0*rdxCp2[1]*bcVals[2]+(12.0*phiUy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-27.0*phiC[0]+72.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.3333333333333333*(3.0*rho[1]*volFac-13.85640646055102*rdxCp2[1]*phiUy[3]-20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiUy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]+25.98076211353316*phiC[0]-207.8460969082653*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiC[3]-60.0*rdxCp2[1]*phiUy[2]+27.0*rdxCp2[0]*phiUx[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiUy[3]-60.0*rdxCp2[0]*phiUx[3]+((-120.0*rdxCp2[1])-300.0*rdxCp2[0])*phiC[3]+51.96152422706631*rdxCp2[0]*phiUx[2]+51.96152422706631*rdxCp2[0]*phiC[2]+(34.64101615137754*phiUy[1]-34.64101615137754*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiC[2]+72.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-27.0*phiC[0])*rdxCp2[1]-13.85640646055102*rdxCp2[0]*phiUx[1]-20.78460969082652*rdxCp2[0]*phiC[1]+(12.0*phiUx[0]-12.0*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiUy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.3333333333333333*(3.0*rho[2]*volFac-13.85640646055102*rdxCp2[0]*phiUx[3]-20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+12.0*rdxCp2[0]*phiUx[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]-207.8460969082653*rdxCp2[1]*bcVals[2]+(25.98076211353316*phiUy[0]+25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[1]*phiUy[3]-40.0*rdxCp2[0]*phiUx[3]+((-300.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+34.64101615137754*rdxCp2[0]*phiUx[2]-34.64101615137754*rdxCp2[0]*phiC[2]+(51.96152422706631*phiUy[1]+51.96152422706631*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac-6.928203230275509*rdxCp2[1]*phiUy[2]-10.39230484541326*rdxCp2[1]*phiC[2]-4.0*rdxCp2[1]*bcVals[2]+(6.0*phiUy[0]-6.0*phiC[0])*rdxCp2[1]-6.928203230275509*rdxCp2[0]*phiUx[1]-10.39230484541326*rdxCp2[0]*phiC[1]+(6.0*phiUx[0]-6.0*phiC[0]-4.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac-41.56921938165305*rdxCp2[1]*phiUy[3]-62.35382907247956*rdxCp2[1]*phiC[3]+(36.0*phiUy[1]-36.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac-41.56921938165305*rdxCp2[0]*phiUx[3]-62.35382907247956*rdxCp2[0]*phiC[3]-60.0*rdxCp2[1]*phiUy[2]+36.0*rdxCp2[0]*phiUx[2]+((-180.0*rdxCp2[1])-36.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiUy[3]-40.0*rdxCp2[0]*phiUx[3]+((-120.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+34.64101615137754*rdxCp2[0]*phiUx[2]-34.64101615137754*rdxCp2[0]*phiC[2]+(34.64101615137754*phiUy[1]-34.64101615137754*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac+24.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]-3.464101615137754*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]-3.464101615137754*rdxCp2[0]*phiUx[1]+3.464101615137754*rdxCp2[0]*phiC[1]+(3.0*phiUx[0]-9.0*phiC[0]+24.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac+3.464101615137754*rdxCp2[1]*phiLy[3]-3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiUx[1]-50.0*rdxCp2[0]*phiC[1]+(8.660254037844386*phiUx[0]+8.660254037844386*phiC[0]-69.28203230275508*bcVals[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac-3.464101615137754*rdxCp2[0]*phiUx[3]+3.464101615137754*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+3.0*rdxCp2[0]*phiUx[2]-10.0*rdxCp2[1]*phiLy[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+((-8.660254037844386*phiLy[0])-8.660254037844386*phiC[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[0]*phiUx[3]-20.0*rdxCp2[1]*phiLy[3]+((-100.0*rdxCp2[1])-100.0*rdxCp2[0])*phiC[3]+17.32050807568877*rdxCp2[0]*phiUx[2]+17.32050807568877*rdxCp2[0]*phiC[2]+((-17.32050807568877*phiLy[1])-17.32050807568877*phiC[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac+8.0*rdxCp2[1]*bcVals[3]+13.85640646055102*rdxCp2[1]*phiLy[2]+20.78460969082652*rdxCp2[1]*phiC[2]+(12.0*phiLy[0]-12.0*phiC[0])*rdxCp2[1]-10.39230484541326*rdxCp2[0]*phiUx[1]+10.39230484541326*rdxCp2[0]*phiC[1]+(9.0*phiUx[0]-27.0*phiC[0]+72.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.3333333333333333*(3.0*rho[1]*volFac+13.85640646055102*rdxCp2[1]*phiLy[3]+20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiLy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiUx[1]-150.0*rdxCp2[0]*phiC[1]+(25.98076211353316*phiUx[0]+25.98076211353316*phiC[0]-207.8460969082653*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac-31.17691453623978*rdxCp2[0]*phiUx[3]+31.17691453623978*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+27.0*rdxCp2[0]*phiUx[2]-60.0*rdxCp2[1]*phiLy[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[0]*phiUx[3]-40.0*rdxCp2[1]*phiLy[3]+((-120.0*rdxCp2[1])-300.0*rdxCp2[0])*phiC[3]+51.96152422706631*rdxCp2[0]*phiUx[2]+51.96152422706631*rdxCp2[0]*phiC[2]+(34.64101615137754*phiC[1]-34.64101615137754*phiLy[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac+72.0*rdxCp2[1]*bcVals[3]+10.39230484541326*rdxCp2[1]*phiLy[2]-10.39230484541326*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-27.0*phiC[0])*rdxCp2[1]-13.85640646055102*rdxCp2[0]*phiUx[1]-20.78460969082652*rdxCp2[0]*phiC[1]+(12.0*phiUx[0]-12.0*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac+31.17691453623978*rdxCp2[1]*phiLy[3]-31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiLy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.3333333333333333*(3.0*rho[2]*volFac-13.85640646055102*rdxCp2[0]*phiUx[3]-20.78460969082652*rdxCp2[0]*phiC[3]+207.8460969082653*rdxCp2[1]*bcVals[3]+12.0*rdxCp2[0]*phiUx[2]-30.0*rdxCp2[1]*phiLy[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+((-25.98076211353316*phiLy[0])-25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[0]*phiUx[3]-60.0*rdxCp2[1]*phiLy[3]+((-300.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+34.64101615137754*rdxCp2[0]*phiUx[2]-34.64101615137754*rdxCp2[0]*phiC[2]+((-51.96152422706631*phiLy[1])-51.96152422706631*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac+4.0*rdxCp2[1]*bcVals[3]+6.928203230275509*rdxCp2[1]*phiLy[2]+10.39230484541326*rdxCp2[1]*phiC[2]+(6.0*phiLy[0]-6.0*phiC[0])*rdxCp2[1]-6.928203230275509*rdxCp2[0]*phiUx[1]-10.39230484541326*rdxCp2[0]*phiC[1]+(6.0*phiUx[0]-6.0*phiC[0]-4.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac+41.56921938165305*rdxCp2[1]*phiLy[3]+62.35382907247956*rdxCp2[1]*phiC[3]+(36.0*phiLy[1]-36.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiUx[1]-180.0*rdxCp2[0]*phiC[1]+(51.96152422706631*phiUx[0]-51.96152422706631*phiC[0]+69.28203230275508*bcVals[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac-41.56921938165305*rdxCp2[0]*phiUx[3]-62.35382907247956*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]+36.0*rdxCp2[0]*phiUx[2]-60.0*rdxCp2[1]*phiLy[2]+((-180.0*rdxCp2[1])-36.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[0]*phiUx[3]-40.0*rdxCp2[1]*phiLy[3]+((-120.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]+34.64101615137754*rdxCp2[0]*phiUx[2]-34.64101615137754*rdxCp2[0]*phiC[2]+(34.64101615137754*phiC[1]-34.64101615137754*phiLy[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac-3.464101615137754*rdxCp2[1]*phiUy[2]+3.464101615137754*rdxCp2[1]*phiC[2]+24.0*rdxCp2[1]*bcVals[2]+(3.0*phiUy[0]-9.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]-3.464101615137754*rdxCp2[0]*phiC[1]+24.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac-3.464101615137754*rdxCp2[1]*phiUy[3]+3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiUy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiLx[1]-50.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+((-8.660254037844386*phiLx[0])-8.660254037844386*phiC[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac+3.464101615137754*rdxCp2[0]*phiLx[3]-3.464101615137754*rdxCp2[0]*phiC[3]-10.0*rdxCp2[1]*phiUy[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]-69.28203230275508*rdxCp2[1]*bcVals[2]+(8.660254037844386*phiUy[0]+8.660254037844386*phiC[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[1]*phiUy[3]-20.0*rdxCp2[0]*phiLx[3]+((-100.0*rdxCp2[1])-100.0*rdxCp2[0])*phiC[3]-17.32050807568877*rdxCp2[0]*phiLx[2]-17.32050807568877*rdxCp2[0]*phiC[2]+(17.32050807568877*phiUy[1]+17.32050807568877*phiC[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-13.85640646055102*rdxCp2[1]*phiUy[2]-20.78460969082652*rdxCp2[1]*phiC[2]-8.0*rdxCp2[1]*bcVals[2]+(12.0*phiUy[0]-12.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiLx[1]-10.39230484541326*rdxCp2[0]*phiC[1]+72.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-27.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.3333333333333333*(3.0*rho[1]*volFac-13.85640646055102*rdxCp2[1]*phiUy[3]-20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiUy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiLx[1]-150.0*rdxCp2[0]*phiC[1]+207.8460969082653*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiLx[0])-25.98076211353316*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac+31.17691453623978*rdxCp2[0]*phiLx[3]-31.17691453623978*rdxCp2[0]*phiC[3]-60.0*rdxCp2[1]*phiUy[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiUy[3]-60.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-300.0*rdxCp2[0])*phiC[3]-51.96152422706631*rdxCp2[0]*phiLx[2]-51.96152422706631*rdxCp2[0]*phiC[2]+(34.64101615137754*phiUy[1]-34.64101615137754*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac-10.39230484541326*rdxCp2[1]*phiUy[2]+10.39230484541326*rdxCp2[1]*phiC[2]+72.0*rdxCp2[1]*bcVals[2]+(9.0*phiUy[0]-27.0*phiC[0])*rdxCp2[1]+13.85640646055102*rdxCp2[0]*phiLx[1]+20.78460969082652*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(12.0*phiLx[0]-12.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac-31.17691453623978*rdxCp2[1]*phiUy[3]+31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiUy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.3333333333333333*(3.0*rho[2]*volFac+13.85640646055102*rdxCp2[0]*phiLx[3]+20.78460969082652*rdxCp2[0]*phiC[3]-30.0*rdxCp2[1]*phiUy[2]+12.0*rdxCp2[0]*phiLx[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]-207.8460969082653*rdxCp2[1]*bcVals[2]+(25.98076211353316*phiUy[0]+25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[1]*phiUy[3]-40.0*rdxCp2[0]*phiLx[3]+((-300.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]-34.64101615137754*rdxCp2[0]*phiLx[2]+34.64101615137754*rdxCp2[0]*phiC[2]+(51.96152422706631*phiUy[1]+51.96152422706631*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac-6.928203230275509*rdxCp2[1]*phiUy[2]-10.39230484541326*rdxCp2[1]*phiC[2]-4.0*rdxCp2[1]*bcVals[2]+(6.0*phiUy[0]-6.0*phiC[0])*rdxCp2[1]+6.928203230275509*rdxCp2[0]*phiLx[1]+10.39230484541326*rdxCp2[0]*phiC[1]+4.0*rdxCp2[0]*bcVals[1]+(6.0*phiLx[0]-6.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac-41.56921938165305*rdxCp2[1]*phiUy[3]-62.35382907247956*rdxCp2[1]*phiC[3]+(36.0*phiUy[1]-36.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac+41.56921938165305*rdxCp2[0]*phiLx[3]+62.35382907247956*rdxCp2[0]*phiC[3]-60.0*rdxCp2[1]*phiUy[2]+36.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-36.0*rdxCp2[0])*phiC[2]+69.28203230275508*rdxCp2[1]*bcVals[2]+(51.96152422706631*phiUy[0]-51.96152422706631*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiUy[3]-40.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]-34.64101615137754*rdxCp2[0]*phiLx[2]+34.64101615137754*rdxCp2[0]*phiC[2]+(34.64101615137754*phiUy[1]-34.64101615137754*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.5*(2.0*rho[0]*volFac+24.0*rdxCp2[1]*bcVals[3]+3.464101615137754*rdxCp2[1]*phiLy[2]-3.464101615137754*rdxCp2[1]*phiC[2]+(3.0*phiLy[0]-9.0*phiC[0])*rdxCp2[1]+3.464101615137754*rdxCp2[0]*phiLx[1]-3.464101615137754*rdxCp2[0]*phiC[1]+24.0*rdxCp2[0]*bcVals[1]+(3.0*phiLx[0]-9.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = rho[1]*volFac+3.464101615137754*rdxCp2[1]*phiLy[3]-3.464101615137754*rdxCp2[1]*phiC[3]+(3.0*phiLy[1]-9.0*phiC[1])*rdxCp2[1]-10.0*rdxCp2[0]*phiLx[1]-50.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+((-8.660254037844386*phiLx[0])-8.660254037844386*phiC[0])*rdxCp2[0]; 
  resOut[2] = rho[2]*volFac+3.464101615137754*rdxCp2[0]*phiLx[3]-3.464101615137754*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]-10.0*rdxCp2[1]*phiLy[2]+3.0*rdxCp2[0]*phiLx[2]+((-50.0*rdxCp2[1])-9.0*rdxCp2[0])*phiC[2]+((-8.660254037844386*phiLy[0])-8.660254037844386*phiC[0])*rdxCp2[1]; 
  resOut[3] = rho[3]*volFac-20.0*rdxCp2[1]*phiLy[3]-20.0*rdxCp2[0]*phiLx[3]+((-100.0*rdxCp2[1])-100.0*rdxCp2[0])*phiC[3]-17.32050807568877*rdxCp2[0]*phiLx[2]-17.32050807568877*rdxCp2[0]*phiC[2]+((-17.32050807568877*phiLy[1])-17.32050807568877*phiC[1])*rdxCp2[1]; 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac+8.0*rdxCp2[1]*bcVals[3]+13.85640646055102*rdxCp2[1]*phiLy[2]+20.78460969082652*rdxCp2[1]*phiC[2]+(12.0*phiLy[0]-12.0*phiC[0])*rdxCp2[1]+10.39230484541326*rdxCp2[0]*phiLx[1]-10.39230484541326*rdxCp2[0]*phiC[1]+72.0*rdxCp2[0]*bcVals[1]+(9.0*phiLx[0]-27.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.3333333333333333*(3.0*rho[1]*volFac+13.85640646055102*rdxCp2[1]*phiLy[3]+20.78460969082652*rdxCp2[1]*phiC[3]+(12.0*phiLy[1]-12.0*phiC[1])*rdxCp2[1]-30.0*rdxCp2[0]*phiLx[1]-150.0*rdxCp2[0]*phiC[1]+207.8460969082653*rdxCp2[0]*bcVals[1]+((-25.98076211353316*phiLx[0])-25.98076211353316*phiC[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac+31.17691453623978*rdxCp2[0]*phiLx[3]-31.17691453623978*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]-60.0*rdxCp2[1]*phiLy[2]+27.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-81.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiLy[3]-60.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-300.0*rdxCp2[0])*phiC[3]-51.96152422706631*rdxCp2[0]*phiLx[2]-51.96152422706631*rdxCp2[0]*phiC[2]+(34.64101615137754*phiC[1]-34.64101615137754*phiLy[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.1666666666666667*(6.0*rho[0]*volFac+72.0*rdxCp2[1]*bcVals[3]+10.39230484541326*rdxCp2[1]*phiLy[2]-10.39230484541326*rdxCp2[1]*phiC[2]+(9.0*phiLy[0]-27.0*phiC[0])*rdxCp2[1]+13.85640646055102*rdxCp2[0]*phiLx[1]+20.78460969082652*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(12.0*phiLx[0]-12.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac+31.17691453623978*rdxCp2[1]*phiLy[3]-31.17691453623978*rdxCp2[1]*phiC[3]+(27.0*phiLy[1]-81.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.3333333333333333*(3.0*rho[2]*volFac+13.85640646055102*rdxCp2[0]*phiLx[3]+20.78460969082652*rdxCp2[0]*phiC[3]+207.8460969082653*rdxCp2[1]*bcVals[3]-30.0*rdxCp2[1]*phiLy[2]+12.0*rdxCp2[0]*phiLx[2]+((-150.0*rdxCp2[1])-12.0*rdxCp2[0])*phiC[2]+((-25.98076211353316*phiLy[0])-25.98076211353316*phiC[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-60.0*rdxCp2[1]*phiLy[3]-40.0*rdxCp2[0]*phiLx[3]+((-300.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]-34.64101615137754*rdxCp2[0]*phiLx[2]+34.64101615137754*rdxCp2[0]*phiC[2]+((-51.96152422706631*phiLy[1])-51.96152422706631*phiC[1])*rdxCp2[1]); 

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
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 
  rdxCp2[1]  = volFac/(dxC[1]*dxC[1]); 
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

  resOut[0] = 0.3333333333333333*(3.0*rho[0]*volFac+4.0*rdxCp2[1]*bcVals[3]+6.928203230275509*rdxCp2[1]*phiLy[2]+10.39230484541326*rdxCp2[1]*phiC[2]+(6.0*phiLy[0]-6.0*phiC[0])*rdxCp2[1]+6.928203230275509*rdxCp2[0]*phiLx[1]+10.39230484541326*rdxCp2[0]*phiC[1]+4.0*rdxCp2[0]*bcVals[1]+(6.0*phiLx[0]-6.0*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1111111111111111*(9.0*rho[1]*volFac+41.56921938165305*rdxCp2[1]*phiLy[3]+62.35382907247956*rdxCp2[1]*phiC[3]+(36.0*phiLy[1]-36.0*phiC[1])*rdxCp2[1]-60.0*rdxCp2[0]*phiLx[1]-180.0*rdxCp2[0]*phiC[1]+69.28203230275508*rdxCp2[0]*bcVals[1]+(51.96152422706631*phiC[0]-51.96152422706631*phiLx[0])*rdxCp2[0]); 
  resOut[2] = 0.1111111111111111*(9.0*rho[2]*volFac+41.56921938165305*rdxCp2[0]*phiLx[3]+62.35382907247956*rdxCp2[0]*phiC[3]+69.28203230275508*rdxCp2[1]*bcVals[3]-60.0*rdxCp2[1]*phiLy[2]+36.0*rdxCp2[0]*phiLx[2]+((-180.0*rdxCp2[1])-36.0*rdxCp2[0])*phiC[2]+(51.96152422706631*phiC[0]-51.96152422706631*phiLy[0])*rdxCp2[1]); 
  resOut[3] = 0.3333333333333333*(3.0*rho[3]*volFac-40.0*rdxCp2[1]*phiLy[3]-40.0*rdxCp2[0]*phiLx[3]+((-120.0*rdxCp2[1])-120.0*rdxCp2[0])*phiC[3]-34.64101615137754*rdxCp2[0]*phiLx[2]+34.64101615137754*rdxCp2[0]*phiC[2]+(34.64101615137754*phiC[1]-34.64101615137754*phiLy[1])*rdxCp2[1]); 

}

