#ifndef EULER_MOD_DELC_H 
#define EULER_MOD_DELC_H 
#include <cmath> 
extern "C" { 
typedef struct {
   double _gasGamma; /* Gas constant */
   double _rpTime; /* Time spent in RP */
   int _numWaves; /* Number of waves in system */
   int _rpType; /* Type of RP to use */
   double _fl[6], _fr[6]; /* Storage for left/right fluxes ([6] as we want to index from 1) */
} EulerEq_t;
} 
#endif 
