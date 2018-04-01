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
void accumulateConfToPhase3x3vMax_P3(const double fact, const double *fconf, double *fphase) 
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
  fphase[28] += 2.828427124746191*fconf[10]*fact; 
  fphase[48] += 2.828427124746191*fconf[11]*fact; 
  fphase[49] += 2.828427124746191*fconf[12]*fact; 
  fphase[50] += 2.828427124746191*fconf[13]*fact; 
  fphase[51] += 2.828427124746191*fconf[14]*fact; 
  fphase[52] += 2.828427124746191*fconf[15]*fact; 
  fphase[53] += 2.828427124746191*fconf[16]*fact; 
  fphase[78] += 2.828427124746191*fconf[17]*fact; 
  fphase[79] += 2.828427124746191*fconf[18]*fact; 
  fphase[80] += 2.828427124746191*fconf[19]*fact; 
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
void assignConfToPhase3x3vMax_P3(const double fact, const double *fconf, double *fphase) 
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
  fphase[28] = 2.828427124746191*fconf[10]*fact; 
  fphase[29] = 0.0; 
  fphase[30] = 0.0; 
  fphase[31] = 0.0; 
  fphase[32] = 0.0; 
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
  fphase[48] = 2.828427124746191*fconf[11]*fact; 
  fphase[49] = 2.828427124746191*fconf[12]*fact; 
  fphase[50] = 2.828427124746191*fconf[13]*fact; 
  fphase[51] = 2.828427124746191*fconf[14]*fact; 
  fphase[52] = 2.828427124746191*fconf[15]*fact; 
  fphase[53] = 2.828427124746191*fconf[16]*fact; 
  fphase[54] = 0.0; 
  fphase[55] = 0.0; 
  fphase[56] = 0.0; 
  fphase[57] = 0.0; 
  fphase[58] = 0.0; 
  fphase[59] = 0.0; 
  fphase[60] = 0.0; 
  fphase[61] = 0.0; 
  fphase[62] = 0.0; 
  fphase[63] = 0.0; 
  fphase[64] = 0.0; 
  fphase[65] = 0.0; 
  fphase[66] = 0.0; 
  fphase[67] = 0.0; 
  fphase[68] = 0.0; 
  fphase[69] = 0.0; 
  fphase[70] = 0.0; 
  fphase[71] = 0.0; 
  fphase[72] = 0.0; 
  fphase[73] = 0.0; 
  fphase[74] = 0.0; 
  fphase[75] = 0.0; 
  fphase[76] = 0.0; 
  fphase[77] = 0.0; 
  fphase[78] = 2.828427124746191*fconf[17]*fact; 
  fphase[79] = 2.828427124746191*fconf[18]*fact; 
  fphase[80] = 2.828427124746191*fconf[19]*fact; 
  fphase[81] = 0.0; 
  fphase[82] = 0.0; 
  fphase[83] = 0.0; 
} 
