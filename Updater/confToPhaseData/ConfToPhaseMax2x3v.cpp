#include <ConfToPhaseModDecl.h>
void accumulateConfToPhase2x3vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[2] += 2.828427124746191*fconf[2]*fact; 
} 
void accumulateConfToPhase2x3vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[2] += 2.828427124746191*fconf[2]*fact; 
  fphase[6] += 2.828427124746191*fconf[3]*fact; 
  fphase[16] += 2.828427124746191*fconf[4]*fact; 
  fphase[17] += 2.828427124746191*fconf[5]*fact; 
} 
void accumulateConfToPhase2x3vMax_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.828427124746191*fconf[0]*fact; 
  fphase[1] += 2.828427124746191*fconf[1]*fact; 
  fphase[2] += 2.828427124746191*fconf[2]*fact; 
  fphase[6] += 2.828427124746191*fconf[3]*fact; 
  fphase[16] += 2.828427124746191*fconf[4]*fact; 
  fphase[17] += 2.828427124746191*fconf[5]*fact; 
  fphase[31] += 2.828427124746191*fconf[6]*fact; 
  fphase[32] += 2.828427124746191*fconf[7]*fact; 
  fphase[51] += 2.828427124746191*fconf[8]*fact; 
  fphase[52] += 2.828427124746191*fconf[9]*fact; 
} 
void assignConfToPhase2x3vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 2.828427124746191*fconf[2]*fact; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
} 
void assignConfToPhase2x3vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 2.828427124746191*fconf[2]*fact; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 2.828427124746191*fconf[3]*fact; 
  fphase[7] = 0.0; 
  fphase[8] = 0.0; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 0.0; 
  fphase[12] = 0.0; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
  fphase[15] = 0.0; 
  fphase[16] = 2.828427124746191*fconf[4]*fact; 
  fphase[17] = 2.828427124746191*fconf[5]*fact; 
  fphase[18] = 0.0; 
  fphase[19] = 0.0; 
  fphase[20] = 0.0; 
} 
void assignConfToPhase2x3vMax_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.828427124746191*fconf[0]*fact; 
  fphase[1] = 2.828427124746191*fconf[1]*fact; 
  fphase[2] = 2.828427124746191*fconf[2]*fact; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 2.828427124746191*fconf[3]*fact; 
  fphase[7] = 0.0; 
  fphase[8] = 0.0; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 0.0; 
  fphase[12] = 0.0; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
  fphase[15] = 0.0; 
  fphase[16] = 2.828427124746191*fconf[4]*fact; 
  fphase[17] = 2.828427124746191*fconf[5]*fact; 
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
  fphase[31] = 2.828427124746191*fconf[6]*fact; 
  fphase[32] = 2.828427124746191*fconf[7]*fact; 
  fphase[33] = 0.0; 
  fphase[34] = 0.0; 
  fphase[35] = 0.0; 
  fphase[36] = 0.0; 
  fphase[37] = 0.0; 
  fphase[38] = 0.0; 
  fphase[39] = 0.0; 
  fphase[40] = 0.0; 
  fphase[41] = 0.0; 
  fphase[42] = 0.0; 
  fphase[43] = 0.0; 
  fphase[44] = 0.0; 
  fphase[45] = 0.0; 
  fphase[46] = 0.0; 
  fphase[47] = 0.0; 
  fphase[48] = 0.0; 
  fphase[49] = 0.0; 
  fphase[50] = 0.0; 
  fphase[51] = 2.828427124746191*fconf[8]*fact; 
  fphase[52] = 2.828427124746191*fconf[9]*fact; 
  fphase[53] = 0.0; 
  fphase[54] = 0.0; 
  fphase[55] = 0.0; 
} 
