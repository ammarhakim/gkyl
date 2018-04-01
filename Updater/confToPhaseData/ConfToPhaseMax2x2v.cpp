#include <ConfToPhaseModDecl.h>
void accumulateConfToPhase2x2vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.0*fconf[0]*fact; 
  fphase[1] += 2.0*fconf[1]*fact; 
  fphase[2] += 2.0*fconf[2]*fact; 
} 
void accumulateConfToPhase2x2vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.0*fconf[0]*fact; 
  fphase[1] += 2.0*fconf[1]*fact; 
  fphase[2] += 2.0*fconf[2]*fact; 
  fphase[5] += 2.0*fconf[3]*fact; 
  fphase[11] += 2.0*fconf[4]*fact; 
  fphase[12] += 2.0*fconf[5]*fact; 
} 
void accumulateConfToPhase2x2vMax_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.0*fconf[0]*fact; 
  fphase[1] += 2.0*fconf[1]*fact; 
  fphase[2] += 2.0*fconf[2]*fact; 
  fphase[5] += 2.0*fconf[3]*fact; 
  fphase[11] += 2.0*fconf[4]*fact; 
  fphase[12] += 2.0*fconf[5]*fact; 
  fphase[19] += 2.0*fconf[6]*fact; 
  fphase[20] += 2.0*fconf[7]*fact; 
  fphase[31] += 2.0*fconf[8]*fact; 
  fphase[32] += 2.0*fconf[9]*fact; 
} 
void assignConfToPhase2x2vMax_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.0*fconf[0]*fact; 
  fphase[1] = 2.0*fconf[1]*fact; 
  fphase[2] = 2.0*fconf[2]*fact; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
} 
void assignConfToPhase2x2vMax_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.0*fconf[0]*fact; 
  fphase[1] = 2.0*fconf[1]*fact; 
  fphase[2] = 2.0*fconf[2]*fact; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 2.0*fconf[3]*fact; 
  fphase[6] = 0.0; 
  fphase[7] = 0.0; 
  fphase[8] = 0.0; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 2.0*fconf[4]*fact; 
  fphase[12] = 2.0*fconf[5]*fact; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
} 
void assignConfToPhase2x2vMax_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.0*fconf[0]*fact; 
  fphase[1] = 2.0*fconf[1]*fact; 
  fphase[2] = 2.0*fconf[2]*fact; 
  fphase[3] = 0.0; 
  fphase[4] = 0.0; 
  fphase[5] = 2.0*fconf[3]*fact; 
  fphase[6] = 0.0; 
  fphase[7] = 0.0; 
  fphase[8] = 0.0; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 2.0*fconf[4]*fact; 
  fphase[12] = 2.0*fconf[5]*fact; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
  fphase[15] = 0.0; 
  fphase[16] = 0.0; 
  fphase[17] = 0.0; 
  fphase[18] = 0.0; 
  fphase[19] = 2.0*fconf[6]*fact; 
  fphase[20] = 2.0*fconf[7]*fact; 
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
  fphase[31] = 2.0*fconf[8]*fact; 
  fphase[32] = 2.0*fconf[9]*fact; 
  fphase[33] = 0.0; 
  fphase[34] = 0.0; 
} 
