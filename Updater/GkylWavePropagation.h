// Gkyl ------------------------------------------------------------------------
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_WAVE_PROPAGATION_H
#define GKYL_WAVE_PROPAGATION_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylCartField.h>
#include <GkylEquationFv.h>

extern "C" 
{
  typedef struct {
      int updateDirs[6];
      int numUpdateDirs;
      double dt;
      double _cfl;
      double _cflm;
      GkylEquationFv_t *equation;
      GkylCartField_t *dtByCell;
      double *buf;
      int sharedMemSize;
  } GkylWavePropagation_t;

  void wavePropagationInitOnDevice(
      const int meqn, const int mwave, const int numBlocks, const int numThreads,
      GkylWavePropagation_t *hyper);

  void wavePropagationDestroyOnDevice(GkylWavePropagation_t *hyper);

  void wavePropagationAdvanceOnDevice(
    const int meqn, const int mwave, const int numBlocks, const int numThreads,
    GkylWavePropagation_t *hyper, GkylCartField_t *qIn, GkylCartField_t *qOut);

  void setDt(GkylWavePropagation_t *hyper, double dt);
}

#endif // GKYL_WAVE_PROPAGATION_H
