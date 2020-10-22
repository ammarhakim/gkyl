#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBC1xMaxP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 8.485281374238571*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-4.898979485566357*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC1xMaxP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 8.485281374238571*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 4.898979485566357*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC1xMaxP2_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 31.1126983722081*fBC-35.77708763999664*fSkin[2]+34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = (-24.49489742783179*fBC)+30.98386676965934*fSkin[2]-29.0*fSkin[1]+17.32050807568877*fSkin[0]; 
  fGhost[2] = 6.324555320336761*fBC-9.0*fSkin[2]+7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC1xMaxP2_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 31.1126983722081*fBC-35.77708763999664*fSkin[2]-34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = 24.49489742783179*fBC-30.98386676965934*fSkin[2]-29.0*fSkin[1]-17.32050807568877*fSkin[0]; 
  fGhost[2] = 6.324555320336761*fBC-9.0*fSkin[2]-7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC1xMaxP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-1.414213562373095*dx*fpBC; 
  fGhost[1] = 0.8164965809277261*dx*fpBC-1.0*fSkin[1]; 

};

void ConstDiffusionBC1xMaxP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[2]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[2]:  DG coefficients in ghost cell.

  fGhost[0] = 1.414213562373095*dx*fpBC+fSkin[0]; 
  fGhost[1] = 0.8164965809277261*dx*fpBC-1.0*fSkin[1]; 

};

void ConstDiffusionBC1xMaxP2_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = (-2.592724864350674*dx*fpBC)-11.18033988749895*fSkin[2]+2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 2.041241452319315*dx*fpBC+11.61895003862225*fSkin[2]-4.0*fSkin[1]; 
  fGhost[2] = (-0.5270462766947298*dx*fpBC)-4.0*fSkin[2]+1.290994448735806*fSkin[1]; 

};

void ConstDiffusionBC1xMaxP2_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 2.592724864350674*dx*fpBC-11.18033988749895*fSkin[2]-2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 2.041241452319315*dx*fpBC-11.61895003862225*fSkin[2]-4.0*fSkin[1]; 
  fGhost[2] = 0.5270462766947298*dx*fpBC-4.0*fSkin[2]-1.290994448735806*fSkin[1]; 

};

void ConstDiffusionBC2xMaxP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-6.928203230275509*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 
  fGhost[2] = -5.0*fSkin[2]; 

};

void ConstDiffusionBC2xMaxP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 6.928203230275509*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 
  fGhost[2] = -5.0*fSkin[2]; 

};

void ConstDiffusionBC2xMaxP1_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC+6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = (-6.928203230275509*fBC)-5.0*fSkin[2]+3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC2xMaxP1_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 12.0*fBC-6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = 6.928203230275509*fBC-5.0*fSkin[2]-3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC2xMaxP2_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = 44.0*fBC-35.77708763999664*fSkin[4]+34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fBC)+30.98386676965934*fSkin[4]-29.0*fSkin[1]+17.32050807568877*fSkin[0]; 
  fGhost[2] = 34.64101615137754*fSkin[3]-21.0*fSkin[2]; 
  fGhost[3] = 17.32050807568877*fSkin[2]-29.0*fSkin[3]; 
  fGhost[4] = 8.94427190999916*fBC-9.0*fSkin[4]+7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 
  fGhost[5] = -21.0*fSkin[5]; 

};

void ConstDiffusionBC2xMaxP2_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = 44.0*fBC-35.77708763999664*fSkin[4]-34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fBC-30.98386676965934*fSkin[4]-29.0*fSkin[1]-17.32050807568877*fSkin[0]; 
  fGhost[2] = (-34.64101615137754*fSkin[3])-21.0*fSkin[2]; 
  fGhost[3] = (-29.0*fSkin[3])-17.32050807568877*fSkin[2]; 
  fGhost[4] = 8.94427190999916*fBC-9.0*fSkin[4]-7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 
  fGhost[5] = -21.0*fSkin[5]; 

};

void ConstDiffusionBC2xMaxP2_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = 44.0*fBC-35.77708763999664*fSkin[5]+34.64101615137754*fSkin[2]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fSkin[3]-21.0*fSkin[1]; 
  fGhost[2] = (-34.64101615137754*fBC)+30.98386676965934*fSkin[5]-29.0*fSkin[2]+17.32050807568877*fSkin[0]; 
  fGhost[3] = 17.32050807568877*fSkin[1]-29.0*fSkin[3]; 
  fGhost[4] = -21.0*fSkin[4]; 
  fGhost[5] = 8.94427190999916*fBC-9.0*fSkin[5]+7.745966692414834*fSkin[2]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC2xMaxP2_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = 44.0*fBC-35.77708763999664*fSkin[5]-34.64101615137754*fSkin[2]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fSkin[3])-21.0*fSkin[1]; 
  fGhost[2] = 34.64101615137754*fBC-30.98386676965934*fSkin[5]-29.0*fSkin[2]-17.32050807568877*fSkin[0]; 
  fGhost[3] = (-29.0*fSkin[3])-17.32050807568877*fSkin[1]; 
  fGhost[4] = -21.0*fSkin[4]; 
  fGhost[5] = 8.94427190999916*fBC-9.0*fSkin[5]-7.745966692414834*fSkin[2]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC2xMaxP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.0*dx*fpBC; 
  fGhost[1] = 1.154700538379252*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 

};

void ConstDiffusionBC2xMaxP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 2.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = 1.154700538379252*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 

};

void ConstDiffusionBC2xMaxP1_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.154700538379252*dx*fpBC-1.0*fSkin[2]; 

};

void ConstDiffusionBC2xMaxP1_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[3]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[3]:  DG coefficients in ghost cell.

  fGhost[0] = 2.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.154700538379252*dx*fpBC-1.0*fSkin[2]; 

};

void ConstDiffusionBC2xMaxP2_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = (-3.666666666666667*dx*fpBC)-11.18033988749895*fSkin[4]+2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*dx*fpBC+11.61895003862225*fSkin[4]-4.0*fSkin[1]; 
  fGhost[2] = 2.886751345948129*fSkin[3]+fSkin[2]; 
  fGhost[3] = -4.0*fSkin[3]; 
  fGhost[4] = (-0.7453559924999299*dx*fpBC)-4.0*fSkin[4]+1.290994448735806*fSkin[1]; 
  fGhost[5] = fSkin[5]; 

};

void ConstDiffusionBC2xMaxP2_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = 3.666666666666667*dx*fpBC-11.18033988749895*fSkin[4]-2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*dx*fpBC-11.61895003862225*fSkin[4]-4.0*fSkin[1]; 
  fGhost[2] = fSkin[2]-2.886751345948129*fSkin[3]; 
  fGhost[3] = -4.0*fSkin[3]; 
  fGhost[4] = 0.7453559924999299*dx*fpBC-4.0*fSkin[4]-1.290994448735806*fSkin[1]; 
  fGhost[5] = fSkin[5]; 

};

void ConstDiffusionBC2xMaxP2_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = (-3.666666666666667*dx*fpBC)-11.18033988749895*fSkin[5]+2.886751345948129*fSkin[2]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*fSkin[3]+fSkin[1]; 
  fGhost[2] = 2.886751345948129*dx*fpBC+11.61895003862225*fSkin[5]-4.0*fSkin[2]; 
  fGhost[3] = -4.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = (-0.7453559924999299*dx*fpBC)-4.0*fSkin[5]+1.290994448735806*fSkin[2]; 

};

void ConstDiffusionBC2xMaxP2_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[6]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[6]:  DG coefficients in ghost cell.

  fGhost[0] = 3.666666666666667*dx*fpBC-11.18033988749895*fSkin[5]-2.886751345948129*fSkin[2]+fSkin[0]; 
  fGhost[1] = fSkin[1]-2.886751345948129*fSkin[3]; 
  fGhost[2] = 2.886751345948129*dx*fpBC-11.61895003862225*fSkin[5]-4.0*fSkin[2]; 
  fGhost[3] = -4.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = 0.7453559924999299*dx*fpBC-4.0*fSkin[5]-1.290994448735806*fSkin[2]; 

};

void ConstDiffusionBC3xMaxP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-9.797958971132715*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = -5.0*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 9.797958971132715*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = -5.0*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC+6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = (-9.797958971132715*fBC)-5.0*fSkin[2]+3.464101615137754*fSkin[0]; 
  fGhost[3] = -5.0*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC-6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = 9.797958971132715*fBC-5.0*fSkin[2]-3.464101615137754*fSkin[0]; 
  fGhost[3] = -5.0*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Dirichlet_X3lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC+6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = (-9.797958971132715*fBC)-5.0*fSkin[3]+3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC3xMaxP1_Dirichlet_X3upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 16.97056274847715*fBC-6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = 9.797958971132715*fBC-5.0*fSkin[3]-3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC3xMaxP2_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 62.22539674441621*fBC-35.77708763999664*fSkin[7]+34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = (-48.98979485566358*fBC)+30.98386676965934*fSkin[7]-29.0*fSkin[1]+17.32050807568877*fSkin[0]; 
  fGhost[2] = 34.64101615137754*fSkin[4]-21.0*fSkin[2]; 
  fGhost[3] = 34.64101615137754*fSkin[5]-21.0*fSkin[3]; 
  fGhost[4] = 17.32050807568877*fSkin[2]-29.0*fSkin[4]; 
  fGhost[5] = 17.32050807568877*fSkin[3]-29.0*fSkin[5]; 
  fGhost[6] = -21.0*fSkin[6]; 
  fGhost[7] = 12.64911064067352*fBC-9.0*fSkin[7]+7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = -21.0*fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 62.22539674441621*fBC-35.77708763999664*fSkin[7]-34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = 48.98979485566358*fBC-30.98386676965934*fSkin[7]-29.0*fSkin[1]-17.32050807568877*fSkin[0]; 
  fGhost[2] = (-34.64101615137754*fSkin[4])-21.0*fSkin[2]; 
  fGhost[3] = (-34.64101615137754*fSkin[5])-21.0*fSkin[3]; 
  fGhost[4] = (-29.0*fSkin[4])-17.32050807568877*fSkin[2]; 
  fGhost[5] = (-29.0*fSkin[5])-17.32050807568877*fSkin[3]; 
  fGhost[6] = -21.0*fSkin[6]; 
  fGhost[7] = 12.64911064067352*fBC-9.0*fSkin[7]-7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = -21.0*fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 62.22539674441621*fBC-35.77708763999664*fSkin[8]+34.64101615137754*fSkin[2]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fSkin[4]-21.0*fSkin[1]; 
  fGhost[2] = (-48.98979485566358*fBC)+30.98386676965934*fSkin[8]-29.0*fSkin[2]+17.32050807568877*fSkin[0]; 
  fGhost[3] = 34.64101615137754*fSkin[6]-21.0*fSkin[3]; 
  fGhost[4] = 17.32050807568877*fSkin[1]-29.0*fSkin[4]; 
  fGhost[5] = -21.0*fSkin[5]; 
  fGhost[6] = 17.32050807568877*fSkin[3]-29.0*fSkin[6]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = 12.64911064067352*fBC-9.0*fSkin[8]+7.745966692414834*fSkin[2]-4.47213595499958*fSkin[0]; 
  fGhost[9] = -21.0*fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 62.22539674441621*fBC-35.77708763999664*fSkin[8]-34.64101615137754*fSkin[2]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fSkin[4])-21.0*fSkin[1]; 
  fGhost[2] = 48.98979485566358*fBC-30.98386676965934*fSkin[8]-29.0*fSkin[2]-17.32050807568877*fSkin[0]; 
  fGhost[3] = (-34.64101615137754*fSkin[6])-21.0*fSkin[3]; 
  fGhost[4] = (-29.0*fSkin[4])-17.32050807568877*fSkin[1]; 
  fGhost[5] = -21.0*fSkin[5]; 
  fGhost[6] = (-29.0*fSkin[6])-17.32050807568877*fSkin[3]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = 12.64911064067352*fBC-9.0*fSkin[8]-7.745966692414834*fSkin[2]-4.47213595499958*fSkin[0]; 
  fGhost[9] = -21.0*fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Dirichlet_X3lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 62.22539674441621*fBC-35.77708763999664*fSkin[9]+34.64101615137754*fSkin[3]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fSkin[5]-21.0*fSkin[1]; 
  fGhost[2] = 34.64101615137754*fSkin[6]-21.0*fSkin[2]; 
  fGhost[3] = (-48.98979485566358*fBC)+30.98386676965934*fSkin[9]-29.0*fSkin[3]+17.32050807568877*fSkin[0]; 
  fGhost[4] = -21.0*fSkin[4]; 
  fGhost[5] = 17.32050807568877*fSkin[1]-29.0*fSkin[5]; 
  fGhost[6] = 17.32050807568877*fSkin[2]-29.0*fSkin[6]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = 12.64911064067352*fBC-9.0*fSkin[9]+7.745966692414834*fSkin[3]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC3xMaxP2_Dirichlet_X3upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 62.22539674441621*fBC-35.77708763999664*fSkin[9]-34.64101615137754*fSkin[3]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fSkin[5])-21.0*fSkin[1]; 
  fGhost[2] = (-34.64101615137754*fSkin[6])-21.0*fSkin[2]; 
  fGhost[3] = 48.98979485566358*fBC-30.98386676965934*fSkin[9]-29.0*fSkin[3]-17.32050807568877*fSkin[0]; 
  fGhost[4] = -21.0*fSkin[4]; 
  fGhost[5] = (-29.0*fSkin[5])-17.32050807568877*fSkin[1]; 
  fGhost[6] = (-29.0*fSkin[6])-17.32050807568877*fSkin[2]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = 12.64911064067352*fBC-9.0*fSkin[9]-7.745966692414834*fSkin[3]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC3xMaxP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.828427124746191*dx*fpBC; 
  fGhost[1] = 1.632993161855453*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*dx*fpBC+fSkin[0]; 
  fGhost[1] = 1.632993161855453*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.828427124746191*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.632993161855453*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 1.632993161855453*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Neumann_X3lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-2.828427124746191*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 1.632993161855453*dx*fpBC-1.0*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP1_Neumann_X3upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[4]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[4]:  DG coefficients in ghost cell.

  fGhost[0] = 2.828427124746191*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 1.632993161855453*dx*fpBC-1.0*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP2_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = (-5.185449728701348*dx*fpBC)-11.18033988749895*fSkin[7]+2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 4.082482904638631*dx*fpBC+11.61895003862225*fSkin[7]-4.0*fSkin[1]; 
  fGhost[2] = 2.886751345948129*fSkin[4]+fSkin[2]; 
  fGhost[3] = 2.886751345948129*fSkin[5]+fSkin[3]; 
  fGhost[4] = -4.0*fSkin[4]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = (-1.05409255338946*dx*fpBC)-4.0*fSkin[7]+1.290994448735806*fSkin[1]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 5.185449728701348*dx*fpBC-11.18033988749895*fSkin[7]-2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 4.082482904638631*dx*fpBC-11.61895003862225*fSkin[7]-4.0*fSkin[1]; 
  fGhost[2] = fSkin[2]-2.886751345948129*fSkin[4]; 
  fGhost[3] = fSkin[3]-2.886751345948129*fSkin[5]; 
  fGhost[4] = -4.0*fSkin[4]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = 1.05409255338946*dx*fpBC-4.0*fSkin[7]-1.290994448735806*fSkin[1]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = (-5.185449728701348*dx*fpBC)-11.18033988749895*fSkin[8]+2.886751345948129*fSkin[2]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*fSkin[4]+fSkin[1]; 
  fGhost[2] = 4.082482904638631*dx*fpBC+11.61895003862225*fSkin[8]-4.0*fSkin[2]; 
  fGhost[3] = 2.886751345948129*fSkin[6]+fSkin[3]; 
  fGhost[4] = -4.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = (-1.05409255338946*dx*fpBC)-4.0*fSkin[8]+1.290994448735806*fSkin[2]; 
  fGhost[9] = fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 5.185449728701348*dx*fpBC-11.18033988749895*fSkin[8]-2.886751345948129*fSkin[2]+fSkin[0]; 
  fGhost[1] = fSkin[1]-2.886751345948129*fSkin[4]; 
  fGhost[2] = 4.082482904638631*dx*fpBC-11.61895003862225*fSkin[8]-4.0*fSkin[2]; 
  fGhost[3] = fSkin[3]-2.886751345948129*fSkin[6]; 
  fGhost[4] = -4.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = 1.05409255338946*dx*fpBC-4.0*fSkin[8]-1.290994448735806*fSkin[2]; 
  fGhost[9] = fSkin[9]; 

};

void ConstDiffusionBC3xMaxP2_Neumann_X3lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = (-5.185449728701348*dx*fpBC)-11.18033988749895*fSkin[9]+2.886751345948129*fSkin[3]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*fSkin[5]+fSkin[1]; 
  fGhost[2] = 2.886751345948129*fSkin[6]+fSkin[2]; 
  fGhost[3] = 4.082482904638631*dx*fpBC+11.61895003862225*fSkin[9]-4.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = (-1.05409255338946*dx*fpBC)-4.0*fSkin[9]+1.290994448735806*fSkin[3]; 

};

void ConstDiffusionBC3xMaxP2_Neumann_X3upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[10]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[10]:  DG coefficients in ghost cell.

  fGhost[0] = 5.185449728701348*dx*fpBC-11.18033988749895*fSkin[9]-2.886751345948129*fSkin[3]+fSkin[0]; 
  fGhost[1] = fSkin[1]-2.886751345948129*fSkin[5]; 
  fGhost[2] = fSkin[2]-2.886751345948129*fSkin[6]; 
  fGhost[3] = 4.082482904638631*dx*fpBC-11.61895003862225*fSkin[9]-4.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = 1.05409255338946*dx*fpBC-4.0*fSkin[9]-1.290994448735806*fSkin[3]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = (-13.85640646055102*fBC)-5.0*fSkin[1]+3.464101615137754*fSkin[0]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = -5.0*fSkin[3]; 
  fGhost[4] = -5.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[1]-5.0*fSkin[0]; 
  fGhost[1] = 13.85640646055102*fBC-5.0*fSkin[1]-3.464101615137754*fSkin[0]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = -5.0*fSkin[3]; 
  fGhost[4] = -5.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = (-13.85640646055102*fBC)-5.0*fSkin[2]+3.464101615137754*fSkin[0]; 
  fGhost[3] = -5.0*fSkin[3]; 
  fGhost[4] = -5.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[2]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = 13.85640646055102*fBC-5.0*fSkin[2]-3.464101615137754*fSkin[0]; 
  fGhost[3] = -5.0*fSkin[3]; 
  fGhost[4] = -5.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X3lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = (-13.85640646055102*fBC)-5.0*fSkin[3]+3.464101615137754*fSkin[0]; 
  fGhost[4] = -5.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X3upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[3]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = 13.85640646055102*fBC-5.0*fSkin[3]-3.464101615137754*fSkin[0]; 
  fGhost[4] = -5.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X4lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC+6.928203230275509*fSkin[4]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = -5.0*fSkin[3]; 
  fGhost[4] = (-13.85640646055102*fBC)-5.0*fSkin[4]+3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC4xMaxP1_Dirichlet_X4upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 24.0*fBC-6.928203230275509*fSkin[4]-5.0*fSkin[0]; 
  fGhost[1] = -5.0*fSkin[1]; 
  fGhost[2] = -5.0*fSkin[2]; 
  fGhost[3] = -5.0*fSkin[3]; 
  fGhost[4] = 13.85640646055102*fBC-5.0*fSkin[4]-3.464101615137754*fSkin[0]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X1lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[11]+34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = (-69.28203230275508*fBC)+30.98386676965934*fSkin[11]-29.0*fSkin[1]+17.32050807568877*fSkin[0]; 
  fGhost[2] = 34.64101615137754*fSkin[5]-21.0*fSkin[2]; 
  fGhost[3] = 34.64101615137754*fSkin[6]-21.0*fSkin[3]; 
  fGhost[4] = 34.64101615137754*fSkin[8]-21.0*fSkin[4]; 
  fGhost[5] = 17.32050807568877*fSkin[2]-29.0*fSkin[5]; 
  fGhost[6] = 17.32050807568877*fSkin[3]-29.0*fSkin[6]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = 17.32050807568877*fSkin[4]-29.0*fSkin[8]; 
  fGhost[9] = -21.0*fSkin[9]; 
  fGhost[10] = -21.0*fSkin[10]; 
  fGhost[11] = 17.88854381999832*fBC-9.0*fSkin[11]+7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 
  fGhost[12] = -21.0*fSkin[12]; 
  fGhost[13] = -21.0*fSkin[13]; 
  fGhost[14] = -21.0*fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X1upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[11]-34.64101615137754*fSkin[1]-21.0*fSkin[0]; 
  fGhost[1] = 69.28203230275508*fBC-30.98386676965934*fSkin[11]-29.0*fSkin[1]-17.32050807568877*fSkin[0]; 
  fGhost[2] = (-34.64101615137754*fSkin[5])-21.0*fSkin[2]; 
  fGhost[3] = (-34.64101615137754*fSkin[6])-21.0*fSkin[3]; 
  fGhost[4] = (-34.64101615137754*fSkin[8])-21.0*fSkin[4]; 
  fGhost[5] = (-29.0*fSkin[5])-17.32050807568877*fSkin[2]; 
  fGhost[6] = (-29.0*fSkin[6])-17.32050807568877*fSkin[3]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = (-29.0*fSkin[8])-17.32050807568877*fSkin[4]; 
  fGhost[9] = -21.0*fSkin[9]; 
  fGhost[10] = -21.0*fSkin[10]; 
  fGhost[11] = 17.88854381999832*fBC-9.0*fSkin[11]-7.745966692414834*fSkin[1]-4.47213595499958*fSkin[0]; 
  fGhost[12] = -21.0*fSkin[12]; 
  fGhost[13] = -21.0*fSkin[13]; 
  fGhost[14] = -21.0*fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X2lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[12]+34.64101615137754*fSkin[2]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fSkin[5]-21.0*fSkin[1]; 
  fGhost[2] = (-69.28203230275508*fBC)+30.98386676965934*fSkin[12]-29.0*fSkin[2]+17.32050807568877*fSkin[0]; 
  fGhost[3] = 34.64101615137754*fSkin[7]-21.0*fSkin[3]; 
  fGhost[4] = 34.64101615137754*fSkin[9]-21.0*fSkin[4]; 
  fGhost[5] = 17.32050807568877*fSkin[1]-29.0*fSkin[5]; 
  fGhost[6] = -21.0*fSkin[6]; 
  fGhost[7] = 17.32050807568877*fSkin[3]-29.0*fSkin[7]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = 17.32050807568877*fSkin[4]-29.0*fSkin[9]; 
  fGhost[10] = -21.0*fSkin[10]; 
  fGhost[11] = -21.0*fSkin[11]; 
  fGhost[12] = 17.88854381999832*fBC-9.0*fSkin[12]+7.745966692414834*fSkin[2]-4.47213595499958*fSkin[0]; 
  fGhost[13] = -21.0*fSkin[13]; 
  fGhost[14] = -21.0*fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X2upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[12]-34.64101615137754*fSkin[2]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fSkin[5])-21.0*fSkin[1]; 
  fGhost[2] = 69.28203230275508*fBC-30.98386676965934*fSkin[12]-29.0*fSkin[2]-17.32050807568877*fSkin[0]; 
  fGhost[3] = (-34.64101615137754*fSkin[7])-21.0*fSkin[3]; 
  fGhost[4] = (-34.64101615137754*fSkin[9])-21.0*fSkin[4]; 
  fGhost[5] = (-29.0*fSkin[5])-17.32050807568877*fSkin[1]; 
  fGhost[6] = -21.0*fSkin[6]; 
  fGhost[7] = (-29.0*fSkin[7])-17.32050807568877*fSkin[3]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = (-29.0*fSkin[9])-17.32050807568877*fSkin[4]; 
  fGhost[10] = -21.0*fSkin[10]; 
  fGhost[11] = -21.0*fSkin[11]; 
  fGhost[12] = 17.88854381999832*fBC-9.0*fSkin[12]-7.745966692414834*fSkin[2]-4.47213595499958*fSkin[0]; 
  fGhost[13] = -21.0*fSkin[13]; 
  fGhost[14] = -21.0*fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X3lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[13]+34.64101615137754*fSkin[3]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fSkin[6]-21.0*fSkin[1]; 
  fGhost[2] = 34.64101615137754*fSkin[7]-21.0*fSkin[2]; 
  fGhost[3] = (-69.28203230275508*fBC)+30.98386676965934*fSkin[13]-29.0*fSkin[3]+17.32050807568877*fSkin[0]; 
  fGhost[4] = 34.64101615137754*fSkin[10]-21.0*fSkin[4]; 
  fGhost[5] = -21.0*fSkin[5]; 
  fGhost[6] = 17.32050807568877*fSkin[1]-29.0*fSkin[6]; 
  fGhost[7] = 17.32050807568877*fSkin[2]-29.0*fSkin[7]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = -21.0*fSkin[9]; 
  fGhost[10] = 17.32050807568877*fSkin[4]-29.0*fSkin[10]; 
  fGhost[11] = -21.0*fSkin[11]; 
  fGhost[12] = -21.0*fSkin[12]; 
  fGhost[13] = 17.88854381999832*fBC-9.0*fSkin[13]+7.745966692414834*fSkin[3]-4.47213595499958*fSkin[0]; 
  fGhost[14] = -21.0*fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X3upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[13]-34.64101615137754*fSkin[3]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fSkin[6])-21.0*fSkin[1]; 
  fGhost[2] = (-34.64101615137754*fSkin[7])-21.0*fSkin[2]; 
  fGhost[3] = 69.28203230275508*fBC-30.98386676965934*fSkin[13]-29.0*fSkin[3]-17.32050807568877*fSkin[0]; 
  fGhost[4] = (-34.64101615137754*fSkin[10])-21.0*fSkin[4]; 
  fGhost[5] = -21.0*fSkin[5]; 
  fGhost[6] = (-29.0*fSkin[6])-17.32050807568877*fSkin[1]; 
  fGhost[7] = (-29.0*fSkin[7])-17.32050807568877*fSkin[2]; 
  fGhost[8] = -21.0*fSkin[8]; 
  fGhost[9] = -21.0*fSkin[9]; 
  fGhost[10] = (-29.0*fSkin[10])-17.32050807568877*fSkin[4]; 
  fGhost[11] = -21.0*fSkin[11]; 
  fGhost[12] = -21.0*fSkin[12]; 
  fGhost[13] = 17.88854381999832*fBC-9.0*fSkin[13]-7.745966692414834*fSkin[3]-4.47213595499958*fSkin[0]; 
  fGhost[14] = -21.0*fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X4lower(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[14]+34.64101615137754*fSkin[4]-21.0*fSkin[0]; 
  fGhost[1] = 34.64101615137754*fSkin[8]-21.0*fSkin[1]; 
  fGhost[2] = 34.64101615137754*fSkin[9]-21.0*fSkin[2]; 
  fGhost[3] = 34.64101615137754*fSkin[10]-21.0*fSkin[3]; 
  fGhost[4] = (-69.28203230275508*fBC)+30.98386676965934*fSkin[14]-29.0*fSkin[4]+17.32050807568877*fSkin[0]; 
  fGhost[5] = -21.0*fSkin[5]; 
  fGhost[6] = -21.0*fSkin[6]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = 17.32050807568877*fSkin[1]-29.0*fSkin[8]; 
  fGhost[9] = 17.32050807568877*fSkin[2]-29.0*fSkin[9]; 
  fGhost[10] = 17.32050807568877*fSkin[3]-29.0*fSkin[10]; 
  fGhost[11] = -21.0*fSkin[11]; 
  fGhost[12] = -21.0*fSkin[12]; 
  fGhost[13] = -21.0*fSkin[13]; 
  fGhost[14] = 17.88854381999832*fBC-9.0*fSkin[14]+7.745966692414834*fSkin[4]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC4xMaxP2_Dirichlet_X4upper(const double dx, const double *fSkin, const double fBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fBC:        Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 88.0*fBC-35.77708763999664*fSkin[14]-34.64101615137754*fSkin[4]-21.0*fSkin[0]; 
  fGhost[1] = (-34.64101615137754*fSkin[8])-21.0*fSkin[1]; 
  fGhost[2] = (-34.64101615137754*fSkin[9])-21.0*fSkin[2]; 
  fGhost[3] = (-34.64101615137754*fSkin[10])-21.0*fSkin[3]; 
  fGhost[4] = 69.28203230275508*fBC-30.98386676965934*fSkin[14]-29.0*fSkin[4]-17.32050807568877*fSkin[0]; 
  fGhost[5] = -21.0*fSkin[5]; 
  fGhost[6] = -21.0*fSkin[6]; 
  fGhost[7] = -21.0*fSkin[7]; 
  fGhost[8] = (-29.0*fSkin[8])-17.32050807568877*fSkin[1]; 
  fGhost[9] = (-29.0*fSkin[9])-17.32050807568877*fSkin[2]; 
  fGhost[10] = (-29.0*fSkin[10])-17.32050807568877*fSkin[3]; 
  fGhost[11] = -21.0*fSkin[11]; 
  fGhost[12] = -21.0*fSkin[12]; 
  fGhost[13] = -21.0*fSkin[13]; 
  fGhost[14] = 17.88854381999832*fBC-9.0*fSkin[14]-7.745966692414834*fSkin[4]-4.47213595499958*fSkin[0]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = 2.309401076758503*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = 2.309401076758503*dx*fpBC-1.0*fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 2.309401076758503*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = 2.309401076758503*dx*fpBC-1.0*fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X3lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 2.309401076758503*dx*fpBC-1.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X3upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = 2.309401076758503*dx*fpBC-1.0*fSkin[3]; 
  fGhost[4] = fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X4lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = fSkin[0]-4.0*dx*fpBC; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = 2.309401076758503*dx*fpBC-1.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP1_Neumann_X4upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[5]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[5]:  DG coefficients in ghost cell.

  fGhost[0] = 4.0*dx*fpBC+fSkin[0]; 
  fGhost[1] = fSkin[1]; 
  fGhost[2] = fSkin[2]; 
  fGhost[3] = fSkin[3]; 
  fGhost[4] = 2.309401076758503*dx*fpBC-1.0*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X1lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = (-7.333333333333333*dx*fpBC)-11.18033988749895*fSkin[11]+2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 5.773502691896258*dx*fpBC+11.61895003862225*fSkin[11]-4.0*fSkin[1]; 
  fGhost[2] = 2.886751345948129*fSkin[5]+fSkin[2]; 
  fGhost[3] = 2.886751345948129*fSkin[6]+fSkin[3]; 
  fGhost[4] = 2.886751345948129*fSkin[8]+fSkin[4]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -4.0*fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = (-1.49071198499986*dx*fpBC)-4.0*fSkin[11]+1.290994448735806*fSkin[1]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X1upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 7.333333333333333*dx*fpBC-11.18033988749895*fSkin[11]-2.886751345948129*fSkin[1]+fSkin[0]; 
  fGhost[1] = 5.773502691896258*dx*fpBC-11.61895003862225*fSkin[11]-4.0*fSkin[1]; 
  fGhost[2] = fSkin[2]-2.886751345948129*fSkin[5]; 
  fGhost[3] = fSkin[3]-2.886751345948129*fSkin[6]; 
  fGhost[4] = fSkin[4]-2.886751345948129*fSkin[8]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -4.0*fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = 1.49071198499986*dx*fpBC-4.0*fSkin[11]-1.290994448735806*fSkin[1]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X2lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = (-7.333333333333333*dx*fpBC)-11.18033988749895*fSkin[12]+2.886751345948129*fSkin[2]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*fSkin[5]+fSkin[1]; 
  fGhost[2] = 5.773502691896258*dx*fpBC+11.61895003862225*fSkin[12]-4.0*fSkin[2]; 
  fGhost[3] = 2.886751345948129*fSkin[7]+fSkin[3]; 
  fGhost[4] = 2.886751345948129*fSkin[9]+fSkin[4]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = -4.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = -4.0*fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = (-1.49071198499986*dx*fpBC)-4.0*fSkin[12]+1.290994448735806*fSkin[2]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X2upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 7.333333333333333*dx*fpBC-11.18033988749895*fSkin[12]-2.886751345948129*fSkin[2]+fSkin[0]; 
  fGhost[1] = fSkin[1]-2.886751345948129*fSkin[5]; 
  fGhost[2] = 5.773502691896258*dx*fpBC-11.61895003862225*fSkin[12]-4.0*fSkin[2]; 
  fGhost[3] = fSkin[3]-2.886751345948129*fSkin[7]; 
  fGhost[4] = fSkin[4]-2.886751345948129*fSkin[9]; 
  fGhost[5] = -4.0*fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = -4.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = -4.0*fSkin[9]; 
  fGhost[10] = fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = 1.49071198499986*dx*fpBC-4.0*fSkin[12]-1.290994448735806*fSkin[2]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X3lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = (-7.333333333333333*dx*fpBC)-11.18033988749895*fSkin[13]+2.886751345948129*fSkin[3]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*fSkin[6]+fSkin[1]; 
  fGhost[2] = 2.886751345948129*fSkin[7]+fSkin[2]; 
  fGhost[3] = 5.773502691896258*dx*fpBC+11.61895003862225*fSkin[13]-4.0*fSkin[3]; 
  fGhost[4] = 2.886751345948129*fSkin[10]+fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = -4.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = -4.0*fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = (-1.49071198499986*dx*fpBC)-4.0*fSkin[13]+1.290994448735806*fSkin[3]; 
  fGhost[14] = fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X3upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 7.333333333333333*dx*fpBC-11.18033988749895*fSkin[13]-2.886751345948129*fSkin[3]+fSkin[0]; 
  fGhost[1] = fSkin[1]-2.886751345948129*fSkin[6]; 
  fGhost[2] = fSkin[2]-2.886751345948129*fSkin[7]; 
  fGhost[3] = 5.773502691896258*dx*fpBC-11.61895003862225*fSkin[13]-4.0*fSkin[3]; 
  fGhost[4] = fSkin[4]-2.886751345948129*fSkin[10]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = -4.0*fSkin[6]; 
  fGhost[7] = -4.0*fSkin[7]; 
  fGhost[8] = fSkin[8]; 
  fGhost[9] = fSkin[9]; 
  fGhost[10] = -4.0*fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = 1.49071198499986*dx*fpBC-4.0*fSkin[13]-1.290994448735806*fSkin[3]; 
  fGhost[14] = fSkin[14]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X4lower(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = (-7.333333333333333*dx*fpBC)-11.18033988749895*fSkin[14]+2.886751345948129*fSkin[4]+fSkin[0]; 
  fGhost[1] = 2.886751345948129*fSkin[8]+fSkin[1]; 
  fGhost[2] = 2.886751345948129*fSkin[9]+fSkin[2]; 
  fGhost[3] = 2.886751345948129*fSkin[10]+fSkin[3]; 
  fGhost[4] = 5.773502691896258*dx*fpBC+11.61895003862225*fSkin[14]-4.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -4.0*fSkin[8]; 
  fGhost[9] = -4.0*fSkin[9]; 
  fGhost[10] = -4.0*fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = (-1.49071198499986*dx*fpBC)-4.0*fSkin[14]+1.290994448735806*fSkin[4]; 

};

void ConstDiffusionBC4xMaxP2_Neumann_X4upper(const double dx, const double *fSkin, const double fpBC, double *fGhost) 
{ 
  // dx:         Cell spacing.
  // fSkin[15]:   DG coefficients in skin cell.
  // fpBC:       Dirichlet boundary value.
  // fGhost[15]:  DG coefficients in ghost cell.

  fGhost[0] = 7.333333333333333*dx*fpBC-11.18033988749895*fSkin[14]-2.886751345948129*fSkin[4]+fSkin[0]; 
  fGhost[1] = fSkin[1]-2.886751345948129*fSkin[8]; 
  fGhost[2] = fSkin[2]-2.886751345948129*fSkin[9]; 
  fGhost[3] = fSkin[3]-2.886751345948129*fSkin[10]; 
  fGhost[4] = 5.773502691896258*dx*fpBC-11.61895003862225*fSkin[14]-4.0*fSkin[4]; 
  fGhost[5] = fSkin[5]; 
  fGhost[6] = fSkin[6]; 
  fGhost[7] = fSkin[7]; 
  fGhost[8] = -4.0*fSkin[8]; 
  fGhost[9] = -4.0*fSkin[9]; 
  fGhost[10] = -4.0*fSkin[10]; 
  fGhost[11] = fSkin[11]; 
  fGhost[12] = fSkin[12]; 
  fGhost[13] = fSkin[13]; 
  fGhost[14] = 1.49071198499986*dx*fpBC-4.0*fSkin[14]-1.290994448735806*fSkin[4]; 

};

