#include <MGpoissonModDecl.h> 
 
void MGpoissonESenergyDG1xSer_P1(const double *dx, double *phi, double *out) 
{ 
  // dx:     cell lengths.
  // phiFld: electrostatic potential.
  // out:    energy out.

  const double volFac = 0.5*dx[0]; 

  out[0] += (12.0*(phi[1]*phi[1])*volFac)/(dx[0]*dx[0]); 
}

void MGpoissonESenergyDG1xSer_P2(const double *dx, double *phi, double *out) 
{ 
  // dx:     cell lengths.
  // phiFld: electrostatic potential.
  // out:    energy out.

  const double volFac = 0.5*dx[0]; 

  out[0] += (60.0*(phi[2]*phi[2])*volFac)/(dx[0]*dx[0])+(12.0*(phi[1]*phi[1])*volFac)/(dx[0]*dx[0]); 
}

