// Gkyl ------------------------------------------------------------------------
//
// Wrappers for calling into kernels that compute moments of the distribution
// function.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylCudaConfig.h>
#include <GkylBasisTypes.h>

// Provides a templated wrapper around the low level C-style kernels
// so they can be called more systematically from C++ code.
template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
class MomentModDecl {
  public:
    // Zeroth moment of Vlasov species distribution function.
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out);
    // First moment of Vlasov species distribution function.
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out);
    // Second moment of Vlasov species distribution function.
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out);
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out);
    // Third moment of Vlasov species distribution function.
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out);
};

template <>
class MomentModDecl<1,1,1,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,1,2,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,1,3,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vTensor_M3i_P3(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,2,1,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,2,2,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,2,3,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vTensor_M3i_P3(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,2,1,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,2,2,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,2,3,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vTensor_M3i_P3(w, dxv, f, out);
    }
    
};


template <>
class MomentModDecl<1,3,1,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,3,2,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,3,3,Gkyl::G_TENSOR_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vTensor_M3i_P3(w, dxv, f, out);
    }
    
};


template <>
class MomentModDecl<1,1,1,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,1,2,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,1,3,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x1vSer_M3i_P3(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,2,1,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,2,2,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,2,3,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x2vSer_M3i_P3(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,2,1,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,2,2,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,2,3,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x2vSer_M3i_P3(w, dxv, f, out);
    }
    
};


template <>
class MomentModDecl<1,3,1,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,3,2,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M3i_P2(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<1,3,3,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M0_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M1i_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M2_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M2ij_P3(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc1x3vSer_M3i_P3(w, dxv, f, out);
    }
    
};


template <>
class MomentModDecl<2,3,1,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M3i_P1(w, dxv, f, out);
    }
    
};

template <>
class MomentModDecl<2,3,2,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M0_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M1i_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M2_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M2ij_P2(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc2x3vSer_M3i_P2(w, dxv, f, out);
    }
    
};


template <>
class MomentModDecl<3,3,1,Gkyl::G_SERENDIPITY_C> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc3x3vSer_M0_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc3x3vSer_M1i_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc3x3vSer_M2_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc3x3vSer_M2ij_P1(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc3x3vSer_M3i_P1(w, dxv, f, out);
    }
    
};



