#include <MGpoissonModDecl.h> 
 
void MGpoissonESenergyDG2xSer_P1(const double *dx, double *phi, double *out) 
{ 
  // dx:     cell lengths.
  // phiFld: electrostatic potential.
  // out:    energy out.

  const double volFac = 0.25*dx[0]*dx[1]; 

  out[0] += (12.0*(phi[3]*phi[3])*volFac)/(dx[0]*dx[0])+(12.0*(phi[1]*phi[1])*volFac)/(dx[0]*dx[0]); 
  out[1] += (12.0*(phi[3]*phi[3])*volFac)/(dx[1]*dx[1])+(12.0*(phi[2]*phi[2])*volFac)/(dx[1]*dx[1]); 
}

void MGpoissonESenergyDG2xSer_P2(const double *dx, double *phi, double *out) 
{ 
  // dx:     cell lengths.
  // phiFld: electrostatic potential.
  // out:    energy out.

  const double volFac = 0.25*dx[0]*dx[1]; 

  out[0] += (12.0*(phi[7]*phi[7])*volFac)/(dx[0]*dx[0])+(60.0*(phi[6]*phi[6])*volFac)/(dx[0]*dx[0])+(60.0*(phi[4]*phi[4])*volFac)/(dx[0]*dx[0])+(12.0*(phi[3]*phi[3])*volFac)/(dx[0]*dx[0])+(12.0*(phi[1]*phi[1])*volFac)/(dx[0]*dx[0]); 
  out[1] += (60.0*(phi[7]*phi[7])*volFac)/(dx[1]*dx[1])+(12.0*(phi[6]*phi[6])*volFac)/(dx[1]*dx[1])+(60.0*(phi[5]*phi[5])*volFac)/(dx[1]*dx[1])+(12.0*(phi[3]*phi[3])*volFac)/(dx[1]*dx[1])+(12.0*(phi[2]*phi[2])*volFac)/(dx[1]*dx[1]); 
}

