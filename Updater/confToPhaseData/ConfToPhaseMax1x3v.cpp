#include <ConfToPhaseModDecl.h>
void accumulateConfToPhase1x3vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
} 
void accumulateConfToPhase1x3vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[11] += 2.828427124746191*fconf[2]*fact; 
} 
void accumulateConfToPhase1x3vMax_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[11] += 2.828427124746191*fconf[2]*fact; 
  fphase[31] += 2.828427124746191*fconf[3]*fact; 
} 
void assignConfToPhase1x3vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 0.0; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
} 
void assignConfToPhase1x3vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 0.0; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 0.0; 
  fphase[7] = 0.0; 
  fphase[8] = 0.0; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 2.828427124746191*fconf[2]*fact; 
  fphase[12] = 0.0; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
} 
void assignConfToPhase1x3vMax_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 0.0; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 0.0; 
  fphase[7] = 0.0; 
  fphase[8] = 0.0; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 2.828427124746191*fconf[2]*fact; 
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
  fphase[22] = 0.0; 
  fphase[23] = 0.0; 
  fphase[24] = 0.0; 
  fphase[25] = 0.0; 
  fphase[26] = 0.0; 
  fphase[27] = 0.0; 
  fphase[28] = 0.0; 
  fphase[29] = 0.0; 
  fphase[30] = 0.0; 
  fphase[31] = 2.828427124746191*fconf[3]*fact; 
  fphase[32] = 0.0; 
  fphase[33] = 0.0; 
  fphase[34] = 0.0; 
} 
