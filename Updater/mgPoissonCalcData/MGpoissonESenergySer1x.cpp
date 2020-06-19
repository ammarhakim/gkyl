#include <MGpoissonModDecl.h> 
 
void MGpoissonESenergyDG1xSer_P1(const double *dx, double *phi, double *out) 
{ 
  // dx:     cell lengths.
  // phiFld: electrostatic potential.
  // out:    energy out.

  out[0] += (12.0*(phi[1]*phi[1]))/(dx[0]*dx[0]); 
}

