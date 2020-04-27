#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBC1xSerDirichlet_Xlower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*fBC-1.0*fSkin[0]; 
  fGhost[1] = fSkin[1]; 

};

void ConstDiffusionBC1xSerDirichlet_Xupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*fBC-1.0*fSkin[0]; 
  fGhost[1] = fSkin[1]; 

};

void ConstDiffusionBC1xSerNeumann_Xlower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-1.414213562373095*dx*fpBC; 
  fGhost[1] = 0.8164965809277261*dx*fpBC-1.0*fSkin[1]; 

};

void ConstDiffusionBC1xSerNeumann_Xupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 1.414213562373095*dx*fpBC+fSkin[0]; 
  fGhost[1] = 0.8164965809277261*dx*fpBC-1.0*fSkin[1]; 

};

void ConstDiffusionBC2xSerDirichlet_Xlower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 8.0*fBC-1.0*fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = -1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC2xSerDirichlet_Xupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 8.0*fBC-1.0*fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = -1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC2xSerDirichlet_Ylower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 8.0*fBC-1.0*fSkin[0]; 
  fGhost[1] = -1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC2xSerDirichlet_Yupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 8.0*fBC-1.0*fSkin[0]; 
  fGhost[1] = -1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC2xSerNeumann_Xlower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = (-4.0*dx*fpBC)+11.54700538379252*fSkin[1]-9.0*fSkin[0]; 
  fGhost[1] = 2.309401076758503*dx*fpBC-13.0*fSkin[1]+10.39230484541326*fSkin[0]; 
  fGhost[2] = 11.54700538379252*fSkin[3]-9.0*fSkin[2]; 
  fGhost[3] = 10.39230484541326*fSkin[2]-13.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerNeumann_Xupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC-11.54700538379252*fSkin[1]-9.0*fSkin[0]; 
  fGhost[1] = 2.309401076758503*dx*fpBC-13.0*fSkin[1]-10.39230484541326*fSkin[0]; 
  fGhost[2] = (-11.54700538379252*fSkin[3])-9.0*fSkin[2]; 
  fGhost[3] = (-13.0*fSkin[3])-10.39230484541326*fSkin[2]; 

};

void ConstDiffusionBC2xSerNeumann_Ylower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = (-4.0*dx*fpBC)+11.54700538379252*fSkin[2]-9.0*fSkin[0]; 
  fGhost[1] = 11.54700538379252*fSkin[3]-9.0*fSkin[1]; 
  fGhost[2] = 2.309401076758503*dx*fpBC-13.0*fSkin[2]+10.39230484541326*fSkin[0]; 
  fGhost[3] = 10.39230484541326*fSkin[1]-13.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerNeumann_Yupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC-11.54700538379252*fSkin[2]-9.0*fSkin[0]; 
  fGhost[1] = (-11.54700538379252*fSkin[3])-9.0*fSkin[1]; 
  fGhost[2] = 2.309401076758503*dx*fpBC-13.0*fSkin[2]-10.39230484541326*fSkin[0]; 
  fGhost[3] = (-13.0*fSkin[3])-10.39230484541326*fSkin[1]; 

};

