#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol5xSerP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[2] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs1345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[3] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs12345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[5]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[3] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[4] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3]+rdxFnu[4])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs2345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[3] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[3] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs35(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[1] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs135(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs1235(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[3] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs235(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs45(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[1] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs145(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[2] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs1245(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[3] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs245(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 
  rdxFnu[2] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs5(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs15(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs125(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs25(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[4]/(dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol5xSerP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[2] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs1345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[3] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs12345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[5]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[3] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[4] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3]+rdxFnu[4])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs2345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[3] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[3] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs35(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs135(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs1235(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[3] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs235(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs45(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[1] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs145(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[2] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs1245(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[3] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs245(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 
  rdxFnu[2] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs5(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs15(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs125(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs25(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[4]/(dx[4]*dx[4]*dx[4]*dx[4]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol5xSerP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[5]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 0.5303300858899105*rdxF[0]*f[1]*nu[71]+0.5303300858899105*f[0]*rdxF[0]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[71]+0.5303300858899105*rdxF[0]*f[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[0]*f[6]*nu[71]+0.5303300858899105*rdxF[0]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[0]*f[9]*nu[71]+0.5303300858899105*rdxF[0]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[0]*f[12]*nu[71]+0.5303300858899105*rdxF[0]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[0]*f[2]*nu[71]+0.5303300858899105*rdxF[0]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[0]*f[4]*nu[71]+0.5303300858899105*rdxF[0]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[0]*f[17]*nu[71]+0.5303300858899105*rdxF[0]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[0]*f[5]*nu[71]+0.5303300858899105*rdxF[0]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[0]*f[20]*nu[71]+0.5303300858899105*rdxF[0]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[0]*f[23]*nu[71]+0.5303300858899105*rdxF[0]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[0]*f[10]*nu[71]+0.5303300858899105*rdxF[0]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[0]*f[13]*nu[71]+0.5303300858899105*rdxF[0]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[0]*f[15]*nu[71]+0.5303300858899105*rdxF[0]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[0]*f[28]*nu[71]+0.5303300858899105*rdxF[0]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[0]*f[24]*nu[71]+0.5303300858899105*rdxF[0]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[64])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[3] += 0.5303300858899105*rdxF[0]*f[1]*nu[71]+0.5303300858899105*f[0]*rdxF[0]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[71]+0.5303300858899105*rdxF[0]*f[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[0]*f[6]*nu[71]+0.5303300858899105*rdxF[0]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[0]*f[9]*nu[71]+0.5303300858899105*rdxF[0]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[0]*f[12]*nu[71]+0.5303300858899105*rdxF[0]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[0]*f[2]*nu[71]+0.5303300858899105*rdxF[0]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[0]*f[4]*nu[71]+0.5303300858899105*rdxF[0]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[0]*f[17]*nu[71]+0.5303300858899105*rdxF[0]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[0]*f[5]*nu[71]+0.5303300858899105*rdxF[0]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[0]*f[20]*nu[71]+0.5303300858899105*rdxF[0]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[0]*f[23]*nu[71]+0.5303300858899105*rdxF[0]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[0]*f[10]*nu[71]+0.5303300858899105*rdxF[0]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[0]*f[13]*nu[71]+0.5303300858899105*rdxF[0]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[0]*f[15]*nu[71]+0.5303300858899105*rdxF[0]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[0]*f[28]*nu[71]+0.5303300858899105*rdxF[0]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[0]*f[24]*nu[71]+0.5303300858899105*rdxF[0]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[64]+rdxF[1]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 
  rdxF[2] = 4.0/(dx[4]*dx[4]); 

  out[3] += 0.5303300858899105*rdxF[0]*f[1]*nu[71]+0.5303300858899105*f[0]*rdxF[0]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[71]+0.5303300858899105*rdxF[0]*f[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[0]*f[6]*nu[71]+0.5303300858899105*rdxF[0]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[0]*f[9]*nu[71]+0.5303300858899105*rdxF[0]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[0]*f[12]*nu[71]+0.5303300858899105*rdxF[0]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[0]*f[2]*nu[71]+0.5303300858899105*rdxF[0]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[0]*f[4]*nu[71]+0.5303300858899105*rdxF[0]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[0]*f[17]*nu[71]+0.5303300858899105*rdxF[0]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[0]*f[5]*nu[71]+0.5303300858899105*rdxF[0]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[0]*f[20]*nu[71]+0.5303300858899105*rdxF[0]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[0]*f[23]*nu[71]+0.5303300858899105*rdxF[0]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[0]*f[10]*nu[71]+0.5303300858899105*rdxF[0]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[0]*f[13]*nu[71]+0.5303300858899105*rdxF[0]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[0]*f[15]*nu[71]+0.5303300858899105*rdxF[0]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[0]*f[28]*nu[71]+0.5303300858899105*rdxF[0]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[0]*f[24]*nu[71]+0.5303300858899105*rdxF[0]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[64]+rdxF[1]*0.1767766952966368*nu[96]+rdxF[2]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs1345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 
  rdxF[3] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[64]+rdxF[2]*0.1767766952966368*nu[96]+rdxF[3]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs12345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[5]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 
  rdxF[3] = 4.0/(dx[3]*dx[3]); 
  rdxF[4] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[2]*nu[71]+0.5303300858899105*f[0]*rdxF[2]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[2]*nu[71]+0.5303300858899105*f[1]*rdxF[2]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[2]*f[6]*nu[71]+0.5303300858899105*f[2]*rdxF[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[2]*f[9]*nu[71]+0.5303300858899105*rdxF[2]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[2]*f[12]*nu[71]+0.5303300858899105*rdxF[2]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*f[2]*rdxF[2]*nu[71]+0.5303300858899105*rdxF[2]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[2]*f[4]*nu[71]+0.5303300858899105*rdxF[2]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[2]*f[17]*nu[71]+0.5303300858899105*rdxF[2]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[2]*f[5]*nu[71]+0.5303300858899105*rdxF[2]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[2]*f[20]*nu[71]+0.5303300858899105*rdxF[2]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[2]*f[23]*nu[71]+0.5303300858899105*rdxF[2]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[2]*f[10]*nu[71]+0.5303300858899105*rdxF[2]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[2]*f[13]*nu[71]+0.5303300858899105*rdxF[2]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[2]*f[15]*nu[71]+0.5303300858899105*rdxF[2]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[2]*f[28]*nu[71]+0.5303300858899105*rdxF[2]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[2]*f[24]*nu[71]+0.5303300858899105*rdxF[2]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[64]+rdxF[3]*0.1767766952966368*nu[96]+rdxF[4]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs2345(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 
  rdxF[3] = 4.0/(dx[4]*dx[4]); 

  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[64]+rdxF[2]*0.1767766952966368*nu[96]+rdxF[3]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[64]+rdxF[2]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 
  rdxF[3] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[2]*nu[71]+0.5303300858899105*f[0]*rdxF[2]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[2]*nu[71]+0.5303300858899105*f[1]*rdxF[2]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[2]*f[6]*nu[71]+0.5303300858899105*f[2]*rdxF[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[2]*f[9]*nu[71]+0.5303300858899105*rdxF[2]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[2]*f[12]*nu[71]+0.5303300858899105*rdxF[2]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*f[2]*rdxF[2]*nu[71]+0.5303300858899105*rdxF[2]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[2]*f[4]*nu[71]+0.5303300858899105*rdxF[2]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[2]*f[17]*nu[71]+0.5303300858899105*rdxF[2]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[2]*f[5]*nu[71]+0.5303300858899105*rdxF[2]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[2]*f[20]*nu[71]+0.5303300858899105*rdxF[2]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[2]*f[23]*nu[71]+0.5303300858899105*rdxF[2]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[2]*f[10]*nu[71]+0.5303300858899105*rdxF[2]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[2]*f[13]*nu[71]+0.5303300858899105*rdxF[2]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[2]*f[15]*nu[71]+0.5303300858899105*rdxF[2]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[2]*f[28]*nu[71]+0.5303300858899105*rdxF[2]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[2]*f[24]*nu[71]+0.5303300858899105*rdxF[2]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[64]+rdxF[3]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[64]+rdxF[2]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs35(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 
  rdxF[1] = 4.0/(dx[4]*dx[4]); 

  out[3] += 0.5303300858899105*rdxF[0]*f[1]*nu[71]+0.5303300858899105*f[0]*rdxF[0]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[71]+0.5303300858899105*rdxF[0]*f[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[0]*f[6]*nu[71]+0.5303300858899105*rdxF[0]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[0]*f[9]*nu[71]+0.5303300858899105*rdxF[0]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[0]*f[12]*nu[71]+0.5303300858899105*rdxF[0]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[0]*f[2]*nu[71]+0.5303300858899105*rdxF[0]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[0]*f[4]*nu[71]+0.5303300858899105*rdxF[0]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[0]*f[17]*nu[71]+0.5303300858899105*rdxF[0]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[0]*f[5]*nu[71]+0.5303300858899105*rdxF[0]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[0]*f[20]*nu[71]+0.5303300858899105*rdxF[0]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[0]*f[23]*nu[71]+0.5303300858899105*rdxF[0]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[0]*f[10]*nu[71]+0.5303300858899105*rdxF[0]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[0]*f[13]*nu[71]+0.5303300858899105*rdxF[0]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[0]*f[15]*nu[71]+0.5303300858899105*rdxF[0]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[0]*f[28]*nu[71]+0.5303300858899105*rdxF[0]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[0]*f[24]*nu[71]+0.5303300858899105*rdxF[0]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[64]+rdxF[1]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs135(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[64]+rdxF[2]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs1235(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 
  rdxF[3] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[2]*nu[71]+0.5303300858899105*f[0]*rdxF[2]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[2]*nu[71]+0.5303300858899105*f[1]*rdxF[2]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[2]*f[6]*nu[71]+0.5303300858899105*f[2]*rdxF[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[2]*f[9]*nu[71]+0.5303300858899105*rdxF[2]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[2]*f[12]*nu[71]+0.5303300858899105*rdxF[2]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*f[2]*rdxF[2]*nu[71]+0.5303300858899105*rdxF[2]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[2]*f[4]*nu[71]+0.5303300858899105*rdxF[2]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[2]*f[17]*nu[71]+0.5303300858899105*rdxF[2]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[2]*f[5]*nu[71]+0.5303300858899105*rdxF[2]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[2]*f[20]*nu[71]+0.5303300858899105*rdxF[2]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[2]*f[23]*nu[71]+0.5303300858899105*rdxF[2]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[2]*f[10]*nu[71]+0.5303300858899105*rdxF[2]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[2]*f[13]*nu[71]+0.5303300858899105*rdxF[2]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[2]*f[15]*nu[71]+0.5303300858899105*rdxF[2]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[2]*f[28]*nu[71]+0.5303300858899105*rdxF[2]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[2]*f[24]*nu[71]+0.5303300858899105*rdxF[2]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[64]+rdxF[3]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs235(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[4]*dx[4]); 

  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[64]+rdxF[2]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[64])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.5303300858899105*f[1]*rdxF[2]*nu[71]+0.5303300858899105*f[0]*rdxF[2]*nu[67]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[2]*nu[71]+0.5303300858899105*f[1]*rdxF[2]*nu[67]+0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.5303300858899105*rdxF[2]*f[6]*nu[71]+0.5303300858899105*f[2]*rdxF[2]*nu[67]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.5303300858899105*rdxF[2]*f[9]*nu[71]+0.5303300858899105*rdxF[2]*f[4]*nu[67]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[14] += 0.5303300858899105*rdxF[2]*f[12]*nu[71]+0.5303300858899105*rdxF[2]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*f[2]*rdxF[2]*nu[71]+0.5303300858899105*rdxF[2]*f[6]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[2]*f[4]*nu[71]+0.5303300858899105*rdxF[2]*f[9]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[19] += 0.5303300858899105*rdxF[2]*f[17]*nu[71]+0.5303300858899105*rdxF[2]*f[10]*nu[67]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[2]*f[5]*nu[71]+0.5303300858899105*rdxF[2]*f[12]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[22] += 0.5303300858899105*rdxF[2]*f[20]*nu[71]+0.5303300858899105*rdxF[2]*f[13]*nu[67]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[25] += 0.5303300858899105*rdxF[2]*f[23]*nu[71]+0.5303300858899105*rdxF[2]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[2]*f[10]*nu[71]+0.5303300858899105*rdxF[2]*f[17]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[2]*f[13]*nu[71]+0.5303300858899105*rdxF[2]*f[20]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[2]*f[15]*nu[71]+0.5303300858899105*rdxF[2]*f[23]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[30] += 0.5303300858899105*rdxF[2]*f[28]*nu[71]+0.5303300858899105*rdxF[2]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[2]*f[24]*nu[71]+0.5303300858899105*rdxF[2]*f[28]*nu[67]+0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[64])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[3] += 0.5303300858899105*f[1]*rdxF[1]*nu[71]+0.5303300858899105*f[0]*rdxF[1]*nu[67]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[1]*nu[71]+0.5303300858899105*f[1]*rdxF[1]*nu[67]; 
  out[8] += 0.5303300858899105*rdxF[1]*f[6]*nu[71]+0.5303300858899105*rdxF[1]*f[2]*nu[67]; 
  out[11] += 0.5303300858899105*rdxF[1]*f[9]*nu[71]+0.5303300858899105*rdxF[1]*f[4]*nu[67]; 
  out[14] += 0.5303300858899105*rdxF[1]*f[12]*nu[71]+0.5303300858899105*rdxF[1]*f[5]*nu[67]; 
  out[16] += 0.5303300858899105*rdxF[1]*f[2]*nu[71]+0.5303300858899105*rdxF[1]*f[6]*nu[67]; 
  out[18] += 0.5303300858899105*rdxF[1]*f[4]*nu[71]+0.5303300858899105*rdxF[1]*f[9]*nu[67]; 
  out[19] += 0.5303300858899105*rdxF[1]*f[17]*nu[71]+0.5303300858899105*rdxF[1]*f[10]*nu[67]; 
  out[21] += 0.5303300858899105*rdxF[1]*f[5]*nu[71]+0.5303300858899105*rdxF[1]*f[12]*nu[67]; 
  out[22] += 0.5303300858899105*rdxF[1]*f[20]*nu[71]+0.5303300858899105*rdxF[1]*f[13]*nu[67]; 
  out[25] += 0.5303300858899105*rdxF[1]*f[23]*nu[71]+0.5303300858899105*rdxF[1]*f[15]*nu[67]; 
  out[26] += 0.5303300858899105*rdxF[1]*f[10]*nu[71]+0.5303300858899105*rdxF[1]*f[17]*nu[67]; 
  out[27] += 0.5303300858899105*rdxF[1]*f[13]*nu[71]+0.5303300858899105*rdxF[1]*f[20]*nu[67]; 
  out[29] += 0.5303300858899105*rdxF[1]*f[15]*nu[71]+0.5303300858899105*rdxF[1]*f[23]*nu[67]; 
  out[30] += 0.5303300858899105*rdxF[1]*f[28]*nu[71]+0.5303300858899105*rdxF[1]*f[24]*nu[67]; 
  out[31] += 0.5303300858899105*rdxF[1]*f[24]*nu[71]+0.5303300858899105*rdxF[1]*f[28]*nu[67]; 

  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[64])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[3]*dx[3]); 


  return (rdxF[0]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs45(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[3]*dx[3]); 
  rdxF[1] = 4.0/(dx[4]*dx[4]); 


  return (rdxF[0]*0.1767766952966368*nu[96]+rdxF[1]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs145(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 
  rdxF[2] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[96]+rdxF[2]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs1245(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 
  rdxF[3] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[96]+rdxF[3]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs245(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 
  rdxF[2] = 4.0/(dx[4]*dx[4]); 


  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[96]+rdxF[2]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 


  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[96])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs5(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[4]*dx[4]); 


  return (rdxF[0]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs15(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs125(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[4]*dx[4]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32]+rdxF[2]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs25(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[4]*dx[4]); 


  return (rdxF[0]*0.1767766952966368*nu[32]+rdxF[1]*0.1767766952966368*nu[128])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 0.5303300858899105*rdxF[0]*f[3]*nu[7]+0.5303300858899105*f[0]*rdxF[0]*nu[1]; 
  out[6] += 0.5303300858899105*rdxF[0]*nu[7]*f[8]+0.5303300858899105*rdxF[0]*nu[1]*f[2]; 
  out[7] += 0.5303300858899105*f[0]*rdxF[0]*nu[7]+0.5303300858899105*rdxF[0]*nu[1]*f[3]; 
  out[9] += 0.5303300858899105*rdxF[0]*nu[7]*f[11]+0.5303300858899105*rdxF[0]*nu[1]*f[4]; 
  out[12] += 0.5303300858899105*rdxF[0]*nu[7]*f[14]+0.5303300858899105*rdxF[0]*nu[1]*f[5]; 
  out[16] += 0.5303300858899105*rdxF[0]*nu[1]*f[8]+0.5303300858899105*rdxF[0]*f[2]*nu[7]; 
  out[17] += 0.5303300858899105*rdxF[0]*nu[7]*f[19]+0.5303300858899105*rdxF[0]*nu[1]*f[10]; 
  out[18] += 0.5303300858899105*rdxF[0]*nu[1]*f[11]+0.5303300858899105*rdxF[0]*f[4]*nu[7]; 
  out[20] += 0.5303300858899105*rdxF[0]*nu[7]*f[22]+0.5303300858899105*rdxF[0]*nu[1]*f[13]; 
  out[21] += 0.5303300858899105*rdxF[0]*nu[1]*f[14]+0.5303300858899105*rdxF[0]*f[5]*nu[7]; 
  out[23] += 0.5303300858899105*rdxF[0]*nu[7]*f[25]+0.5303300858899105*rdxF[0]*nu[1]*f[15]; 
  out[26] += 0.5303300858899105*rdxF[0]*nu[1]*f[19]+0.5303300858899105*rdxF[0]*nu[7]*f[10]; 
  out[27] += 0.5303300858899105*rdxF[0]*nu[1]*f[22]+0.5303300858899105*rdxF[0]*nu[7]*f[13]; 
  out[28] += 0.5303300858899105*rdxF[0]*nu[7]*f[30]+0.5303300858899105*rdxF[0]*nu[1]*f[24]; 
  out[29] += 0.5303300858899105*rdxF[0]*nu[1]*f[25]+0.5303300858899105*rdxF[0]*nu[7]*f[15]; 
  out[31] += 0.5303300858899105*rdxF[0]*nu[1]*f[30]+0.5303300858899105*rdxF[0]*nu[7]*f[24]; 

  return (rdxF[0]*0.1767766952966368*nu[0]+rdxF[1]*0.1767766952966368*nu[32])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol5xSerP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates.
  // dx[5]:  Cell spacing.
  // nu[160]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 


  return (rdxF[0]*0.1767766952966368*nu[32])*0.6666666666666666;

} 
