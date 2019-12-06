#ifndef RANGE_DEVICE_IMPL_H
#define RANGE_DEVICE_IMPL_H

#include <cstdlib>

extern "C"
{
    typedef struct {
      int32_t ndim; int32_t lower[6]; int32_t upper[6];
      int rowMajorIndexerCoeff[7], colMajorIndexerCoeff[7];
    } Range_t;
}

namespace Gkyl {
  // In the index() methods below the -1 accounts for the fact that
  // the Range_t object's indexer coefficients assume a Lua 1-based
  // indexing

  /** Provides indexing into a N-dimension box: this class is not
   * templated and perhaps may provide a slightly slower index
   * functions
   */
  class GenIndexer {
    public:
      /** Return linear index given an N-dimensional index into range
       *
       * @param range Range object
       * @param idx NDIM size index into range
       * @return Linear index (0-start)
       */
      __device__ __inline__ int index(const Range_t *range, int *idx) {
        const int *ac = range->rowMajorIndexerCoeff;
        int loc = -1+ac[0];
        for (int i=0; i<range->ndim; ++i)
          loc += idx[i]*ac[i+1];
        return loc;
      }

      /** Computes NDIM index in a range object given a linear index
       * into it (0-based)
       *
       * @param range Range object
       * @param loc Linear location into range
       * @param idx [out] NDIM index object
       */
      __device__ __inline__ void invIndex(const Range_t *range, int loc, int *idx) {
        const int *ac = range->rowMajorIndexerCoeff;
        div_t dout;
        int n = loc;
        for (int i=1; i<=range->ndim; ++i) {
          int quot = n/ac[i];
          int rem = n % ac[i];
          idx[i-1] = quot + range->lower[i-1];
          n = rem;
        }
      }
  };

  /* Class that provides inverse indexer (not for direct use) */
  template <unsigned NDIM>
  class _InvIndexer {
    public:
      /** Computes NDIM index in a range object given a linear index
       * into it (0-based)
       *
       * @param range Range object
       * @param loc Linear location into range
       * @param idx [out] NDIM index object
       */
      __device__ __inline__ void invIndex(const Range_t *range, int loc, int idx[NDIM]) {
        const int *ac = range->rowMajorIndexerCoeff;
        div_t dout;
        int n = loc;
        for (int i=1; i<=NDIM; ++i) {
          int quot = n/ac[i];
          int rem = n % ac[i];
          idx[i-1] = quot + range->lower[i-1];
          n = rem;
        }
      }
  };

  /** Provides indexing into a N-dimension box */
  template <unsigned NDIM>
  class Indexer : public _InvIndexer<NDIM> {
    public:
      __device__  __inline__ int index(const Range_t *range, int idx[NDIM]) {
        const int *ac = range->rowMajorIndexerCoeff;
        int loc = -1+ac[0];
        for (int i=0; i<NDIM; ++i)
          loc += idx[i]*ac[i+1];
        return loc;
      }
  };

  // 1D indexer
  template <>
  class Indexer<1> : public _InvIndexer<1> {
    public:
      __device__ __inline__ int index(const Range_t *range, int i1) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1];
      }
  };

  // 2D indexer
  template <>
  class Indexer<2> {
    public:
      __device__ __inline__ int index(const Range_t *range, int i1, int i2) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2];
      }
  };

  // 3D indexer
  template <>
  class Indexer<3> {
    public:
      __device__ __inline__ int index(const Range_t *range, int i1, int i2, int i3) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3];
      }
  };

  // 4D indexer
  template <>
  class Indexer<4> {
    public:
      __device__ __inline__ int index(const Range_t *range, int i1, int i2, int i3, int i4) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4];
      }
  };

  // 5D indexer
  template <>
  class Indexer<5> {
    public:
      __device__ __inline__ int index(const Range_t *range, int i1, int i2, int i3, int i4, int i5) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5];
      }
  };

  // 6D indexer
  template <>
  class Indexer<6> {
    public:
      __device__ __inline__ int index(const Range_t *range, int i1, int i2, int i3, int i4, int i5, int i6) {
        const int *ac = range->rowMajorIndexerCoeff;
        return -1+ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]+i6*ac[6];
      }
  };
}

#endif // RANGE_DEVICE_IMPL_H
