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

  // In the index() methods below the -1 accounts for the fact that
  // Lua flat layout assumes a starting index of 1

  // 1D indexer
  template <>
  class Indexer<1> {
    public:
      __device__ int index(const Range_t *range, int i1) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1];
      }
  };

  // 2D indexer
  template <>
  class Indexer<2> {
    public:
      __device__ int index(const Range_t *range, int i1, int i2) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2];
      }
  };

  // 3D indexer
  template <>
  class Indexer<3> {
    public:
      __device__ int index(const Range_t *range, int i1, int i2, int i3) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3];
      }
  };

  // 4D indexer
  template <>
  class Indexer<4> {
    public:
      __device__ int index(const Range_t *range, int i1, int i2, int i3, int i4) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4];
      }
  };

  // 5D indexer
  template <>
  class Indexer<5> {
    public:
      __device__ int index(const Range_t *range, int i1, int i2, int i3, int i4, int i5) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5];
      }
  };

  // 6D indexer
  template <>
  class Indexer<6> {
    public:
      __device__ int index(const Range_t *range, int i1, int i2, int i3, int i4, int i5, int i6) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]+i6*ac[6];
      }
  };
}

#endif // RANGE_DEVICE_IMPL_H
