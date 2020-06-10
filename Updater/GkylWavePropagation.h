// Gkyl ------------------------------------------------------------------------
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_WAVE_PROPAGATION_H
#define GKYL_WAVE_PROPAGATION_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylEquation.h>
#include <GkylEuler.h>
#include <GkylCartField.h>

extern "C" 
{
  typedef struct {
      int updateDirs[6];
      int32_t numUpdateDirs;
      double dt;
      Gkyl::Euler *equation;
  } GkylWavePropagation_t;
    
  void wavePropagationAdvanceOnDevice(
      int numBlocks, int numThreads, GkylWavePropagation_t *hyper,
      GkylCartField_t *qIn, GkylCartField_t *qOut);
}

#endif // GKYL_WAVE_PROPAGATION_H
