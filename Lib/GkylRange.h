// Gkyl ------------------------------------------------------------------------
//
// Range and indexer objects
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

// Gkyl includes
#include <GkylCudaConfig.h>

// std includes
#include <cstdlib>
#include <cmath>
#include <stdint.h>

extern "C"
{
    typedef struct {
      int ndim; int lower[6]; int upper[6];
      int rowMajorIndexerCoeff[7], colMajorIndexerCoeff[7];
      __host__ __device__ inline int volume() const {
        int v = 1;
        #pragma unroll
        for (int i=0; i<6; ++i)
          v *= (upper[i]-lower[i]+1);
        return v;
      }
      __host__ __device__ inline int shape(int dir) const {
        return upper[dir]-lower[dir]+1;
      }        
    } GkylRange_t;
}

namespace Gkyl {
  // In the index() methods below the -1 accounts for the fact that
  // the GkylRange_t object's indexer coefficients assume a Lua 1-based
  // indexing

  // Layout specifiers for indexers
  enum class Layout {
    rowMajor = 1,
    colMajor = 2,
  };

  /** Provides indexing into a N-dimension box: this class is not
   * templated and perhaps may provide a slightly slower index
   * functions
   */
  class GenIndexer {
    public:
      /**
       * @param range Range object
       * @param layout Layout of data: Layout::rowMajor, Layout::colMajor
       */
      __host__ __device__ GenIndexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : range(range), layout(layout) {
        ac = (Gkyl::Layout::rowMajor==layout) ? range->rowMajorIndexerCoeff : range->colMajorIndexerCoeff;
      }

      /** Return linear index given an N-dimensional index into range
       *
       * @param idx NDIM size index into range
       * @return Linear index (0-start)
       */
      __host__ __device__ inline int index(const int* __restrict__ idx) {
        int loc = -1+ac[0];
        #pragma unroll
        for (int i=0; i<6; ++i)
          loc += idx[i]*ac[i+1];
        return loc;
      }

      /** Computes NDIM index in a range object given a linear index
       * into it (0-based)
       *
       * @param loc Linear location into range
       * @param idx [out] NDIM index object
       */
      __host__ __device__ inline void invIndex(int loc, int *idx) {
        if (Gkyl::Layout::rowMajor == layout)
          return invIndexRowMajor(loc, idx);
        else
          return invIndexColMajor(loc, idx);
      }

      __host__ __device__ inline int coeff(int dir) {
        return ac[dir];
      }

    private:
      /** Pointer to range object indexed by this indexer */
      const GkylRange_t *range;
      /** Pointer to indexing data */
      const int *ac;
      /** Layout */
      Gkyl::Layout layout;

      /** Row-major inverse indexer */
      __host__ __device__ inline void invIndexRowMajor(int loc, int *idx) {
        int n = loc;
        #pragma unroll
        for (int i=1; i<=6; ++i) {
          int quot = n/ac[i];
          int rem = n % ac[i];
          idx[i-1] = quot + range->lower[i-1];
          n = rem;
        }
      }

      /** Col-major inverse indexer */
      __host__ __device__ inline void invIndexColMajor(int loc, int *idx) {
        int n = loc;
        #pragma unroll
        for (int i=6; i>=1; --i) {
          int quot = n/ac[i];
          int rem = n % ac[i];
          idx[i-1] = quot + range->lower[i-1];
          n = rem;
        }
      }      
  };

  /* Private class that provides NDIM indexer and inverse indexer (not for direct use) */
  template <unsigned NDIM>
  class _Indexer {
    public:
      /** Return linear index given an N-dimensional index into range
       *
       * @param idx NDIM size index into range
       * @return Linear index (0-start)
       */
      __host__ __device__ inline int64_t index(int *idx) {
        int loc = -1+ac[0];
        for (int i=0; i<NDIM; ++i)
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
      __host__ __device__ inline void invIndex(int64_t loc, int idx[NDIM]) {
        if (Gkyl::Layout::rowMajor == layout)
          return invIndexRowMajor(loc, idx);
        else
          return invIndexColMajor(loc, idx);
      }


    protected:
      /** Ctor: Protected so only for use in derived classes 
       * 
       * @param range Range object to index
       */
      __host__ __device__ _Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : range(range), layout(layout) {
        ac = (Gkyl::Layout::rowMajor==layout) ? range->rowMajorIndexerCoeff : range->colMajorIndexerCoeff;
      }

      /** Pointer to range object indexed by this indexer */
      const GkylRange_t *range;
      /** Pointer to indexing data */
      const int *ac;

    private:
      __host__ __device__ inline void invIndexRowMajor(int loc, int idx[NDIM]) {
        int n = loc;
        for (int i=1; i<=NDIM; ++i) {
          int quot = n/ac[i];
          int rem = n % ac[i];
          idx[i-1] = quot + range->lower[i-1];
          n = rem;
        }
      }

      __host__ __device__ inline void invIndexColMajor(int loc, int idx[NDIM]) {
        int n = loc;
        for (int i=NDIM; i>=1; --i) {
          int quot = n/ac[i];
          int rem = n % ac[i];
          idx[i-1] = quot + range->lower[i-1];
          n = rem;
        }
      }

      /** Layout */
      Gkyl::Layout layout;      
  };

  /** Provides indexing into a N-dimension box */
  template <unsigned NDIM>
  class Indexer : public _Indexer<NDIM> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<NDIM>(range, layout) {
      }

      __host__ __device__  inline int64_t index(int idx[NDIM]) {
        int loc = -1+this->ac[0];
        for (int i=0; i<NDIM; ++i)
          loc += idx[i]*this->ac[i+1];
        return loc;
      }
  };

  // 1D indexer
  template <>
  class Indexer<1> : public _Indexer<1> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<1>(range, layout) {
      }

      __host__ __device__ inline int64_t index(int i1) {
        return -1+this->ac[0]+i1*this->ac[1];
      }
  };

  // 2D indexer
  template <>
  class Indexer<2> : public _Indexer<2> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<2>(range, layout) {
      }

      __host__ __device__ inline int64_t index(int i1, int i2) {
        return -1+this->ac[0]+i1*this->ac[1]+i2*this->ac[2];
      }
  };

  // 3D indexer
  template <>
  class Indexer<3> : public _Indexer<3> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<3>(range, layout) {
      }

      __host__ __device__ inline int64_t index(int i1, int i2, int i3) {
        return -1+this->ac[0]+i1*this->ac[1]+i2*this->ac[2]+i3*this->ac[3];
      }
  };

  // 4D indexer
  template <>
  class Indexer<4> : public _Indexer<4> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<4>(range, layout) {
      }

      __host__ __device__ inline int64_t index(int i1, int i2, int i3, int i4) {
        return -1+this->ac[0]+i1*this->ac[1]+i2*this->ac[2]+i3*this->ac[3]+i4*this->ac[4];
      }
  };

  // 5D indexer
  template <>
  class Indexer<5> : public _Indexer<5> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<5>(range, layout) {
      }

      __host__ __device__ inline int64_t index(int i1, int i2, int i3, int i4, int i5) {
        return -1+this->ac[0]+i1*this->ac[1]+i2*this->ac[2]+i3*this->ac[3]+i4*this->ac[4]+i5*this->ac[5];
      }
  };

  // 6D indexer
  template <>
  class Indexer<6> : public _Indexer<6> {
    public:
      __host__ __device__ Indexer(const GkylRange_t *range, Gkyl::Layout layout=Gkyl::Layout::rowMajor)
        : _Indexer<6>(range, layout) {
      }

      __host__ __device__ inline int64_t index(int i1, int i2, int i3, int i4, int i5, int i6) {
        return -1+this->ac[0]+i1*this->ac[1]+i2*this->ac[2]+i3*this->ac[3]+i4*this->ac[4]+i5*this->ac[5]+i6*this->ac[6];
      }
  };

  /**
   * Range object: N-dimensional index set.
   */
  class Range {
    public:
      __host__ __device__ Range(const GkylRange_t *range)
        : range(range) {
      }
      __host__ __device__ inline int ndim() const {
        return range->ndim;
      }
      __host__ __device__ inline int lower(int dir) const {
        return range->lower[dir];
      }
      __host__ __device__ inline int upper(int dir) const {
        return range->upper[dir];
      }
      __host__ __device__ inline int shape(int dir) const {
        return range->shape(dir);
      }
      __host__ __device__ inline int volume() const {
        return range->volume();
      }

      __host__ __device__ inline Gkyl::GenIndexer genIndexer(Gkyl::Layout layout=Gkyl::Layout::rowMajor) {
        return Gkyl::GenIndexer(range, layout);
      }

      template<unsigned NDIM>
      __host__ __device__ inline Gkyl::Indexer<NDIM> indexer(Gkyl::Layout layout=Gkyl::Layout::rowMajor) {
        return Gkyl::Indexer<NDIM>(range, layout);
      }

    private:
      /* Pointer to range object (created elsewhere) */
      const GkylRange_t *range;
  };

  /**
   * Private functions (don't use these directly)
   */
  static void calcRowMajorIndexerCoeff(GkylRange_t& range) {
    int ndim = range.ndim;
    range.rowMajorIndexerCoeff[ndim] = 1;
    for (int i=ndim-1; i>=1; --i)
      range.rowMajorIndexerCoeff[i] = range.rowMajorIndexerCoeff[i+1]*range.shape(i+1-1);
    int start = 0;
    for (int i=0; i<ndim; ++i)
      start += range.rowMajorIndexerCoeff[i+1]*range.lower[i];
    range.rowMajorIndexerCoeff[0] = 1-start;
  }
  
  static void calcColMajorIndexerCoeff(GkylRange_t& range) {
    int ndim = range.ndim;
    range.colMajorIndexerCoeff[1] = 1;
    for (int i=2; i<=ndim; ++i)
      range.colMajorIndexerCoeff[i] = range.colMajorIndexerCoeff[i-1]*range.shape(i-1-1);
    int start = 0;
    for (int i=0; i<ndim; ++i)
      start += range.colMajorIndexerCoeff[i+1]*range.lower[i];
    range.colMajorIndexerCoeff[0] = 1-start;
  }

  /**
   * Function to construct indexer coeff arrays in GkylRange_t
   * struct. This assumes that the ndim, lower and upper are set
   * properly.
   *
   * @param [in/out] On out, row and col indexer coeffs are constructed
   */
  static void calcIndexerCoeff(GkylRange_t& range) {
    calcRowMajorIndexerCoeff(range);
    calcColMajorIndexerCoeff(range);
  }
}
