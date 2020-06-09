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
      double dt;
      Gkyl::Vlasov *equation;
      GkylCartField_t *cflRateByCell;
  } GkylHyperDisCont_t;
    
  void advanceOnDevice(int numBlocks, int numThreads, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut);
  void advanceOnDevice_shared(int numBlocks, int numThreads, int numComponents, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut);
  void setDtAndCflRate(GkylHyperDisCont_t *hyper, double dt, GkylCartField_t *cflRate);
}

void setDtAndCflRate(GkylHyperDisCont_t *hyper, double dt, GkylCartField_t *cflRate) {
  setDtAndCflRateOnDevice<<<1, 1>>>(hyper, emField);
}

__global__ void setDtAndCflRateOnDevice(GkylHyperDisCont_t *hyper, double dt, GkylCartField_t *cflRate) {
  hyper->dt = dt;
  hyper->cflRateByCell = cflRate;
}

#endif // GKYL_HYPER_DISCONT_H
