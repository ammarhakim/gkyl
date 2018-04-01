#include <ConfToPhaseModDecl.h>
void accumulateConfToPhase1x1vSer_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 1.414213562373095*fconf[0]*fact; 
  fphase[1] += 1.414213562373095*fconf[1]*fact; 
} 
void accumulateConfToPhase1x1vSer_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 1.414213562373095*fconf[0]*fact; 
  fphase[1] += 1.414213562373095*fconf[1]*fact; 
  fphase[4] += 1.414213562373095*fconf[2]*fact; 
} 
void accumulateConfToPhase1x1vSer_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 1.414213562373095*fconf[0]*fact; 
  fphase[1] += 1.414213562373095*fconf[1]*fact; 
  fphase[4] += 1.414213562373095*fconf[2]*fact; 
  fphase[8] += 1.414213562373095*fconf[3]*fact; 
} 
void assignConfToPhase1x1vSer_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 1.414213562373095*fconf[0]*fact; 
  fphase[1] = 1.414213562373095*fconf[1]*fact; 
  fphase[2] = 0.0; 
  fphase[3] = 0.0; 
} 
void assignConfToPhase1x1vSer_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 1.414213562373095*fconf[0]*fact; 
  fphase[1] = 1.414213562373095*fconf[1]*fact; 
  fphase[2] = 0.0; 
  fphase[3] = 0.0; 
  fphase[4] = 1.414213562373095*fconf[2]*fact; 
  fphase[5] = 0.0; 
  fphase[6] = 0.0; 
  fphase[7] = 0.0; 
} 
void assignConfToPhase1x1vSer_P3(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 1.414213562373095*fconf[0]*fact; 
  fphase[1] = 1.414213562373095*fconf[1]*fact; 
  fphase[2] = 0.0; 
  fphase[3] = 0.0; 
  fphase[4] = 1.414213562373095*fconf[2]*fact; 
  fphase[5] = 0.0; 
  fphase[6] = 0.0; 
  fphase[7] = 0.0; 
  fphase[8] = 1.414213562373095*fconf[3]*fact; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 0.0; 
} 
