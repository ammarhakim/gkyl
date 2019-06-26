#include <ConfToPhaseModDecl.h>
void accumulateConfToPhase3x3vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[2] += 2.828427124746191*fconf[2]*fact; 
  fphase[3] += 2.828427124746191*fconf[3]*fact; 
} 
void accumulateConfToPhase3x3vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[2] += 2.828427124746191*fconf[2]*fact; 
  fphase[3] += 2.828427124746191*fconf[3]*fact; 
  fphase[7] += 2.828427124746191*fconf[4]*fact; 
  fphase[8] += 2.828427124746191*fconf[5]*fact; 
  fphase[9] += 2.828427124746191*fconf[6]*fact; 
  fphase[22] += 2.828427124746191*fconf[7]*fact; 
  fphase[23] += 2.828427124746191*fconf[8]*fact; 
  fphase[24] += 2.828427124746191*fconf[9]*fact; 
} 
void assignConfToPhase3x3vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 2.828427124746191*fconf[2]*fact; 
  fphase[3] = 2.828427124746191*fconf[3]*fact; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 0.0; 
} 
void assignConfToPhase3x3vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 2.828427124746191*fconf[2]*fact; 
  fphase[3] = 2.828427124746191*fconf[3]*fact; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 0.0; 
  fphase[7] = 2.828427124746191*fconf[4]*fact; 
  fphase[8] = 2.828427124746191*fconf[5]*fact; 
  fphase[9] = 2.828427124746191*fconf[6]*fact; 
  fphase[10] = 0.0; 
  fphase[11] = 0.0; 
  fphase[12] = 0.0; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
  fphase[15] = 0.0; 
  fphase[16] = 0.0; 
  fphase[17] = 0.0; 
  fphase[18] = 0.0; 
  fphase[19] = 0.0; 
  fphase[20] = 0.0; 
  fphase[21] = 0.0; 
  fphase[22] = 2.828427124746191*fconf[7]*fact; 
  fphase[23] = 2.828427124746191*fconf[8]*fact; 
  fphase[24] = 2.828427124746191*fconf[9]*fact; 
  fphase[25] = 0.0; 
  fphase[26] = 0.0; 
  fphase[27] = 0.0; 
}
