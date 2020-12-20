#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBC1xSerP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 8.485281374238571*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-4.898979485566357*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC1xSerP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 8.485281374238571*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 4.898979485566357*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC1xSerP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-1.414213562373095*dx*fpBC; 
  fGhost[1] = 0.8164965809277261*dx*fpBC-1.0*fSkin[1]; 

};

void ConstDiffusionBC1xSerP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 1.414213562373095*dx*fpBC+fSkin[0]; 
  fGhost[1] = 0.8164965809277261*dx*fpBC-1.0*fSkin[1]; 

};

void ConstDiffusionBC2xSerP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 
  fGhost[2] = 6.928203230275509*fSkin[3]-5.0*fSkin[2]; 
  fGhost[3] = 3.464101615137754*fSkin[2]-5.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 
  fGhost[2] = (-6.928203230275509*fSkin[3])-5.0*fSkin[2]; 
  fGhost[3] = (-5.0*fSkin[3])-3.464101615137754*fSkin[2]; 

};

void ConstDiffusionBC2xSerP1_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC+6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fSkin[3]-5.0*fSkin[1]; 
  fGhost[2] = (-6.928203230275509*fBC)-5.0*fSkin[2]+3.464101615137754*fSkin[0]; 
  fGhost[3] = 3.464101615137754*fSkin[1]-5.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerP1_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC-6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fSkin[3])-5.0*fSkin[1]; 
  fGhost[2] = 6.928203230275509*fBC-5.0*fSkin[2]-3.464101615137754*fSkin[0]; 
  fGhost[3] = (-5.0*fSkin[3])-3.464101615137754*fSkin[1]; 

};

void ConstDiffusionBC2xSerP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.0*dx*fpBC; 
  fGhost[1] = 1.154700538379252*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = -1.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 2.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = 1.154700538379252*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = -1.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerP1_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.154700538379252*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = -1.0*fSkin[3]; 

};

void ConstDiffusionBC2xSerP1_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 2.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.154700538379252*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = -1.0*fSkin[3]; 

};

void ConstDiffusionBC3xSerP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-9.797958971132715*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 
  fGhost[2] = 6.928203230275509*fSkin[4]-5.0*fSkin[2]; 
  fGhost[3] = 6.928203230275509*fSkin[5]-5.0*fSkin[3]; 
  fGhost[4] = 3.464101615137754*fSkin[2]-5.0*fSkin[4]; 
  fGhost[5] = 3.464101615137754*fSkin[3]-5.0*fSkin[5]; 
  fGhost[6] = 6.928203230275509*fSkin[7]-5.0*fSkin[6]; 
  fGhost[7] = 3.464101615137754*fSkin[6]-5.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 9.797958971132715*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 
  fGhost[2] = (-6.928203230275509*fSkin[4])-5.0*fSkin[2]; 
  fGhost[3] = (-6.928203230275509*fSkin[5])-5.0*fSkin[3]; 
  fGhost[4] = (-5.0*fSkin[4])-3.464101615137754*fSkin[2]; 
  fGhost[5] = (-5.0*fSkin[5])-3.464101615137754*fSkin[3]; 
  fGhost[6] = (-6.928203230275509*fSkin[7])-5.0*fSkin[6]; 
  fGhost[7] = (-5.0*fSkin[7])-3.464101615137754*fSkin[6]; 

};

void ConstDiffusionBC3xSerP1_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC+6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fSkin[4]-5.0*fSkin[1]; 
  fGhost[2] = (-9.797958971132715*fBC)-5.0*fSkin[2]+3.464101615137754*fSkin[0]; 
  fGhost[3] = 6.928203230275509*fSkin[6]-5.0*fSkin[3]; 
  fGhost[4] = 3.464101615137754*fSkin[1]-5.0*fSkin[4]; 
  fGhost[5] = 6.928203230275509*fSkin[7]-5.0*fSkin[5]; 
  fGhost[6] = 3.464101615137754*fSkin[3]-5.0*fSkin[6]; 
  fGhost[7] = 3.464101615137754*fSkin[5]-5.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC-6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fSkin[4])-5.0*fSkin[1]; 
  fGhost[2] = 9.797958971132715*fBC-5.0*fSkin[2]-3.464101615137754*fSkin[0]; 
  fGhost[3] = (-6.928203230275509*fSkin[6])-5.0*fSkin[3]; 
  fGhost[4] = (-5.0*fSkin[4])-3.464101615137754*fSkin[1]; 
  fGhost[5] = (-6.928203230275509*fSkin[7])-5.0*fSkin[5]; 
  fGhost[6] = (-5.0*fSkin[6])-3.464101615137754*fSkin[3]; 
  fGhost[7] = (-5.0*fSkin[7])-3.464101615137754*fSkin[5]; 

};

void ConstDiffusionBC3xSerP1_Dirichlet_X3lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC+6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fSkin[5]-5.0*fSkin[1]; 
  fGhost[2] = 6.928203230275509*fSkin[6]-5.0*fSkin[2]; 
  fGhost[3] = (-9.797958971132715*fBC)-5.0*fSkin[3]+3.464101615137754*fSkin[0]; 
  fGhost[4] = 6.928203230275509*fSkin[7]-5.0*fSkin[4]; 
  fGhost[5] = 3.464101615137754*fSkin[1]-5.0*fSkin[5]; 
  fGhost[6] = 3.464101615137754*fSkin[2]-5.0*fSkin[6]; 
  fGhost[7] = 3.464101615137754*fSkin[4]-5.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Dirichlet_X3upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC-6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fSkin[5])-5.0*fSkin[1]; 
  fGhost[2] = (-6.928203230275509*fSkin[6])-5.0*fSkin[2]; 
  fGhost[3] = 9.797958971132715*fBC-5.0*fSkin[3]-3.464101615137754*fSkin[0]; 
  fGhost[4] = (-6.928203230275509*fSkin[7])-5.0*fSkin[4]; 
  fGhost[5] = (-5.0*fSkin[5])-3.464101615137754*fSkin[1]; 
  fGhost[6] = (-5.0*fSkin[6])-3.464101615137754*fSkin[2]; 
  fGhost[7] = (-5.0*fSkin[7])-3.464101615137754*fSkin[4]; 

};

void ConstDiffusionBC3xSerP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.828427124746191*dx*fpBC; 
  fGhost[1] = 1.632993161855453*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = -1.0*fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*dx*fpBC+fSkin[0]; 
  fGhost[1] = 1.632993161855453*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = -1.0*fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.828427124746191*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.632993161855453*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = -1.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.632993161855453*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = -1.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Neumann_X3lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.828427124746191*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 1.632993161855453*dx*fpBC-1.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 

};

void ConstDiffusionBC3xSerP1_Neumann_X3upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[8]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[8]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 1.632993161855453*dx*fpBC-1.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-13.85640646055102*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 
  fGhost[2] = 6.928203230275509*fSkin[5]-5.0*fSkin[2]; 
  fGhost[3] = 6.928203230275509*fSkin[6]-5.0*fSkin[3]; 
  fGhost[4] = 6.928203230275509*fSkin[8]-5.0*fSkin[4]; 
  fGhost[5] = 3.464101615137754*fSkin[2]-5.0*fSkin[5]; 
  fGhost[6] = 3.464101615137754*fSkin[3]-5.0*fSkin[6]; 
  fGhost[7] = 6.928203230275509*fSkin[11]-5.0*fSkin[7]; 
  fGhost[8] = 3.464101615137754*fSkin[4]-5.0*fSkin[8]; 
  fGhost[9] = 6.928203230275509*fSkin[12]-5.0*fSkin[9]; 
  fGhost[10] = 6.928203230275509*fSkin[13]-5.0*fSkin[10]; 
  fGhost[11] = 3.464101615137754*fSkin[7]-5.0*fSkin[11]; 
  fGhost[12] = 3.464101615137754*fSkin[9]-5.0*fSkin[12]; 
  fGhost[13] = 3.464101615137754*fSkin[10]-5.0*fSkin[13]; 
  fGhost[14] = 6.928203230275509*fSkin[15]-5.0*fSkin[14]; 
  fGhost[15] = 3.464101615137754*fSkin[14]-5.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 13.85640646055102*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 
  fGhost[2] = (-6.928203230275509*fSkin[5])-5.0*fSkin[2]; 
  fGhost[3] = (-6.928203230275509*fSkin[6])-5.0*fSkin[3]; 
  fGhost[4] = (-6.928203230275509*fSkin[8])-5.0*fSkin[4]; 
  fGhost[5] = (-5.0*fSkin[5])-3.464101615137754*fSkin[2]; 
  fGhost[6] = (-5.0*fSkin[6])-3.464101615137754*fSkin[3]; 
  fGhost[7] = (-6.928203230275509*fSkin[11])-5.0*fSkin[7]; 
  fGhost[8] = (-5.0*fSkin[8])-3.464101615137754*fSkin[4]; 
  fGhost[9] = (-6.928203230275509*fSkin[12])-5.0*fSkin[9]; 
  fGhost[10] = (-6.928203230275509*fSkin[13])-5.0*fSkin[10]; 
  fGhost[11] = (-5.0*fSkin[11])-3.464101615137754*fSkin[7]; 
  fGhost[12] = (-5.0*fSkin[12])-3.464101615137754*fSkin[9]; 
  fGhost[13] = (-5.0*fSkin[13])-3.464101615137754*fSkin[10]; 
  fGhost[14] = (-6.928203230275509*fSkin[15])-5.0*fSkin[14]; 
  fGhost[15] = (-5.0*fSkin[15])-3.464101615137754*fSkin[14]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fSkin[5]-5.0*fSkin[1]; 
  fGhost[2] = (-13.85640646055102*fBC)-5.0*fSkin[2]+3.464101615137754*fSkin[0]; 
  fGhost[3] = 6.928203230275509*fSkin[7]-5.0*fSkin[3]; 
  fGhost[4] = 6.928203230275509*fSkin[9]-5.0*fSkin[4]; 
  fGhost[5] = 3.464101615137754*fSkin[1]-5.0*fSkin[5]; 
  fGhost[6] = 6.928203230275509*fSkin[11]-5.0*fSkin[6]; 
  fGhost[7] = 3.464101615137754*fSkin[3]-5.0*fSkin[7]; 
  fGhost[8] = 6.928203230275509*fSkin[12]-5.0*fSkin[8]; 
  fGhost[9] = 3.464101615137754*fSkin[4]-5.0*fSkin[9]; 
  fGhost[10] = 6.928203230275509*fSkin[14]-5.0*fSkin[10]; 
  fGhost[11] = 3.464101615137754*fSkin[6]-5.0*fSkin[11]; 
  fGhost[12] = 3.464101615137754*fSkin[8]-5.0*fSkin[12]; 
  fGhost[13] = 6.928203230275509*fSkin[15]-5.0*fSkin[13]; 
  fGhost[14] = 3.464101615137754*fSkin[10]-5.0*fSkin[14]; 
  fGhost[15] = 3.464101615137754*fSkin[13]-5.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fSkin[5])-5.0*fSkin[1]; 
  fGhost[2] = 13.85640646055102*fBC-5.0*fSkin[2]-3.464101615137754*fSkin[0]; 
  fGhost[3] = (-6.928203230275509*fSkin[7])-5.0*fSkin[3]; 
  fGhost[4] = (-6.928203230275509*fSkin[9])-5.0*fSkin[4]; 
  fGhost[5] = (-5.0*fSkin[5])-3.464101615137754*fSkin[1]; 
  fGhost[6] = (-6.928203230275509*fSkin[11])-5.0*fSkin[6]; 
  fGhost[7] = (-5.0*fSkin[7])-3.464101615137754*fSkin[3]; 
  fGhost[8] = (-6.928203230275509*fSkin[12])-5.0*fSkin[8]; 
  fGhost[9] = (-5.0*fSkin[9])-3.464101615137754*fSkin[4]; 
  fGhost[10] = (-6.928203230275509*fSkin[14])-5.0*fSkin[10]; 
  fGhost[11] = (-5.0*fSkin[11])-3.464101615137754*fSkin[6]; 
  fGhost[12] = (-5.0*fSkin[12])-3.464101615137754*fSkin[8]; 
  fGhost[13] = (-6.928203230275509*fSkin[15])-5.0*fSkin[13]; 
  fGhost[14] = (-5.0*fSkin[14])-3.464101615137754*fSkin[10]; 
  fGhost[15] = (-5.0*fSkin[15])-3.464101615137754*fSkin[13]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X3lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fSkin[6]-5.0*fSkin[1]; 
  fGhost[2] = 6.928203230275509*fSkin[7]-5.0*fSkin[2]; 
  fGhost[3] = (-13.85640646055102*fBC)-5.0*fSkin[3]+3.464101615137754*fSkin[0]; 
  fGhost[4] = 6.928203230275509*fSkin[10]-5.0*fSkin[4]; 
  fGhost[5] = 6.928203230275509*fSkin[11]-5.0*fSkin[5]; 
  fGhost[6] = 3.464101615137754*fSkin[1]-5.0*fSkin[6]; 
  fGhost[7] = 3.464101615137754*fSkin[2]-5.0*fSkin[7]; 
  fGhost[8] = 6.928203230275509*fSkin[13]-5.0*fSkin[8]; 
  fGhost[9] = 6.928203230275509*fSkin[14]-5.0*fSkin[9]; 
  fGhost[10] = 3.464101615137754*fSkin[4]-5.0*fSkin[10]; 
  fGhost[11] = 3.464101615137754*fSkin[5]-5.0*fSkin[11]; 
  fGhost[12] = 6.928203230275509*fSkin[15]-5.0*fSkin[12]; 
  fGhost[13] = 3.464101615137754*fSkin[8]-5.0*fSkin[13]; 
  fGhost[14] = 3.464101615137754*fSkin[9]-5.0*fSkin[14]; 
  fGhost[15] = 3.464101615137754*fSkin[12]-5.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X3upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fSkin[6])-5.0*fSkin[1]; 
  fGhost[2] = (-6.928203230275509*fSkin[7])-5.0*fSkin[2]; 
  fGhost[3] = 13.85640646055102*fBC-5.0*fSkin[3]-3.464101615137754*fSkin[0]; 
  fGhost[4] = (-6.928203230275509*fSkin[10])-5.0*fSkin[4]; 
  fGhost[5] = (-6.928203230275509*fSkin[11])-5.0*fSkin[5]; 
  fGhost[6] = (-5.0*fSkin[6])-3.464101615137754*fSkin[1]; 
  fGhost[7] = (-5.0*fSkin[7])-3.464101615137754*fSkin[2]; 
  fGhost[8] = (-6.928203230275509*fSkin[13])-5.0*fSkin[8]; 
  fGhost[9] = (-6.928203230275509*fSkin[14])-5.0*fSkin[9]; 
  fGhost[10] = (-5.0*fSkin[10])-3.464101615137754*fSkin[4]; 
  fGhost[11] = (-5.0*fSkin[11])-3.464101615137754*fSkin[5]; 
  fGhost[12] = (-6.928203230275509*fSkin[15])-5.0*fSkin[12]; 
  fGhost[13] = (-5.0*fSkin[13])-3.464101615137754*fSkin[8]; 
  fGhost[14] = (-5.0*fSkin[14])-3.464101615137754*fSkin[9]; 
  fGhost[15] = (-5.0*fSkin[15])-3.464101615137754*fSkin[12]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X4lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[4]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fSkin[8]-5.0*fSkin[1]; 
  fGhost[2] = 6.928203230275509*fSkin[9]-5.0*fSkin[2]; 
  fGhost[3] = 6.928203230275509*fSkin[10]-5.0*fSkin[3]; 
  fGhost[4] = (-13.85640646055102*fBC)-5.0*fSkin[4]+3.464101615137754*fSkin[0]; 
  fGhost[5] = 6.928203230275509*fSkin[12]-5.0*fSkin[5]; 
  fGhost[6] = 6.928203230275509*fSkin[13]-5.0*fSkin[6]; 
  fGhost[7] = 6.928203230275509*fSkin[14]-5.0*fSkin[7]; 
  fGhost[8] = 3.464101615137754*fSkin[1]-5.0*fSkin[8]; 
  fGhost[9] = 3.464101615137754*fSkin[2]-5.0*fSkin[9]; 
  fGhost[10] = 3.464101615137754*fSkin[3]-5.0*fSkin[10]; 
  fGhost[11] = 6.928203230275509*fSkin[15]-5.0*fSkin[11]; 
  fGhost[12] = 3.464101615137754*fSkin[5]-5.0*fSkin[12]; 
  fGhost[13] = 3.464101615137754*fSkin[6]-5.0*fSkin[13]; 
  fGhost[14] = 3.464101615137754*fSkin[7]-5.0*fSkin[14]; 
  fGhost[15] = 3.464101615137754*fSkin[11]-5.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Dirichlet_X4upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[4]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fSkin[8])-5.0*fSkin[1]; 
  fGhost[2] = (-6.928203230275509*fSkin[9])-5.0*fSkin[2]; 
  fGhost[3] = (-6.928203230275509*fSkin[10])-5.0*fSkin[3]; 
  fGhost[4] = 13.85640646055102*fBC-5.0*fSkin[4]-3.464101615137754*fSkin[0]; 
  fGhost[5] = (-6.928203230275509*fSkin[12])-5.0*fSkin[5]; 
  fGhost[6] = (-6.928203230275509*fSkin[13])-5.0*fSkin[6]; 
  fGhost[7] = (-6.928203230275509*fSkin[14])-5.0*fSkin[7]; 
  fGhost[8] = (-5.0*fSkin[8])-3.464101615137754*fSkin[1]; 
  fGhost[9] = (-5.0*fSkin[9])-3.464101615137754*fSkin[2]; 
  fGhost[10] = (-5.0*fSkin[10])-3.464101615137754*fSkin[3]; 
  fGhost[11] = (-6.928203230275509*fSkin[15])-5.0*fSkin[11]; 
  fGhost[12] = (-5.0*fSkin[12])-3.464101615137754*fSkin[5]; 
  fGhost[13] = (-5.0*fSkin[13])-3.464101615137754*fSkin[6]; 
  fGhost[14] = (-5.0*fSkin[14])-3.464101615137754*fSkin[7]; 
  fGhost[15] = (-5.0*fSkin[15])-3.464101615137754*fSkin[11]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = 2.309401076758503*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -1.0*fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = -1.0*fSkin[11]; 
  fGhost[12] = -1.0*fSkin[12]; 
  fGhost[13] = -1.0*fSkin[13]; 
  fGhost[14] = fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = 2.309401076758503*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -1.0*fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = -1.0*fSkin[11]; 
  fGhost[12] = -1.0*fSkin[12]; 
  fGhost[13] = -1.0*fSkin[13]; 
  fGhost[14] = fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 2.309401076758503*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = -1.0*fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = -1.0*fSkin[11]; 
  fGhost[12] = -1.0*fSkin[12]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = -1.0*fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 2.309401076758503*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -1.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = -1.0*fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = -1.0*fSkin[11]; 
  fGhost[12] = -1.0*fSkin[12]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = -1.0*fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X3lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 2.309401076758503*dx*fpBC-1.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = -1.0*fSkin[10]; 
  fGhost[11] = -1.0*fSkin[11]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = -1.0*fSkin[13]; 
  fGhost[14] = -1.0*fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X3upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 2.309401076758503*dx*fpBC-1.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -1.0*fSkin[6]; 
  fGhost[7] = -1.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = -1.0*fSkin[10]; 
  fGhost[11] = -1.0*fSkin[11]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = -1.0*fSkin[13]; 
  fGhost[14] = -1.0*fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X4lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = 2.309401076758503*dx*fpBC-1.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -1.0*fSkin[8]; 
  fGhost[9] = -1.0*fSkin[9]; 
  fGhost[10] = -1.0*fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = -1.0*fSkin[12]; 
  fGhost[13] = -1.0*fSkin[13]; 
  fGhost[14] = -1.0*fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

void ConstDiffusionBC4xSerP1_Neumann_X4upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[16]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[16]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = 2.309401076758503*dx*fpBC-1.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -1.0*fSkin[8]; 
  fGhost[9] = -1.0*fSkin[9]; 
  fGhost[10] = -1.0*fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = -1.0*fSkin[12]; 
  fGhost[13] = -1.0*fSkin[13]; 
  fGhost[14] = -1.0*fSkin[14]; 
  fGhost[15] = -1.0*fSkin[15]; 

};

