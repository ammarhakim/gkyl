#pragma once

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylRange.h>
#include <GkylRectCart.h>
#include <GkylEquation.h>

extern "C" 
{
  typedef struct {
      int updateDirs[6];
      bool zeroFluxFlags[6];
      int32_t numUpdateDirs;
      bool updateVolumeTerm;
      Equation *equation;
      GkylCartField_t *cflRateByCell;
  } GkylHyperDisCont_t;
    
  void advanceOnDevice(int numThreads, int numBlocks, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut);
}

