// Gkyl ------------------------------------------------------------------------
//
// MF 2022/10/09: A dummy header file to make waf work (temporary hack).
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_DUMMY_STRUCT_H
#define GKYL_DUMMY_STRUCT_H

// Gkyl includes
#include <GkylCudaConfig.h>

extern "C"
{
  typedef struct {
      int up;
      bool zeroF;
      double dt;
      double *maxs;
  } GkylDummyStruct_t;

  void advanceOnDevice(const int numBlocks, const int numThreads, const GkylDummyStruct_t *dummy);
  void setDt(GkylDummyStruct_t *dummy);
}

#endif // GKYL_DUMMY_STRUCT_H
