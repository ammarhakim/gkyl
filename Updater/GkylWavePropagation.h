// Gkyl ------------------------------------------------------------------------
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_WAVE_PROPAGATION_H
#define GKYL_WAVE_PROPAGATION_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylEquationFv.h>
#include <GkylCartField.h>

extern "C" 
{
  typedef struct {
      int updateDirs[6];
      int32_t numUpdateDirs;
      double dt;
      double _cfl;
      double _cflm;
      GkylEquationFv_t *equation;
      GkylCartField_t *dtByCell;
  } GkylWavePropagation_t;
    
  void wavePropagationAdvanceOnDevice(
      int numBlocks, int numThreads, GkylWavePropagation_t *hyper,
      GkylCartField_t *qIn, GkylCartField_t *qOut);

  void setDt(GkylWavePropagation_t *hyper, double dt);
}

#endif // GKYL_WAVE_PROPAGATION_H
