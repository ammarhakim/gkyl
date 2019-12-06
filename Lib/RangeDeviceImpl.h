#ifndef RANGE_DEVICE_IMPL_H
#define RANGE_DEVICE_IMPL_H

extern "C" 
{
    typedef struct {
      int32_t ndim; int32_t lower[6]; int32_t upper[6];
      int rowMajorIndexerCoeff[7], colMajorIndexerCoeff[7];
    } Range_t;
}

namespace Gkyl {

  /** Provides indexing into a N-dimension box */
  template <unsigned NDIM>
  class Indexer {
    public:
      __device__ int index(const Range_t *range, int idx[NDIM]);
  };

  // 1D indexer
  template <>
  class Indexer<1> {
    public:
      __device__ int index(const Range_t *range, int i1) {
        return 1;
      }
  };

  // 2D indexer
  template <>
  class Indexer<2> {
    public:
      __device__ int index(const Range_t *range, int i1, int i2) {
        return 2;
      }
  };
}

#endif // RANGE_DEVICE_IMPL_H
