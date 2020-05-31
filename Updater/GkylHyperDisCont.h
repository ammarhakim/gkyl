// Gkyl ------------------------------------------------------------------------
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_HYPER_DISCONT_H
#define GKYL_HYPER_DISCONT_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylEquation.h>
#include <GkylVlasov.h>
#include <GkylCartField.h>

extern "C" 
{
  typedef struct {
      int updateDirs[6];
      bool zeroFluxFlags[6];
      int32_t numUpdateDirs;
      bool updateVolumeTerm;
      Gkyl::Vlasov *equation;
      GkylCartField_t *cflRateByCell;
  } GkylHyperDisCont_t;
    
  void advanceOnDevice(int numThreads, int numBlocks, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut);
}

#endif // GKYL_HYPER_DISCONT_H
