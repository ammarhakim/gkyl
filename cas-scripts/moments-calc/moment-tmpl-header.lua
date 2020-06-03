-- Generate header file for with template for calling moment kernels.

outDir      = '/Users/mana/max-out/'   -- Output directory.
outFileName = 'DistFuncMomentCalcTmpl.h'    -- Specify where to write the output file.

local xsys = require "xsys"

momentHeaderTemplateTopString = [[
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

]]

momentHeaderTemplateBottomString = [[

]]

momentHeaderTemplateString = [[
|for ci = CMIN, CMAX do
| minV = ci
| if VMIN>minV then minV=VMIN end
| for vi = minV, VMAX do
|  for pi = PMIN, PMAX do
template <>
class MomentModDecl<${ci},${vi},${pi},Gkyl::${basisNm}> {
  public:
    __host__ __device__ static void calcM0(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc${ci}x${vi}v${basisShortNm}_M0_P${pi}(w, dxv, f, out);
    }
    __host__ __device__ static void calcM1i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc${ci}x${vi}v${basisShortNm}_M1i_P${pi}(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc${ci}x${vi}v${basisShortNm}_M2_P${pi}(w, dxv, f, out);
    }
    __host__ __device__ static void calcM2ij(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc${ci}x${vi}v${basisShortNm}_M2ij_P${pi}(w, dxv, f, out);
    }
    __host__ __device__ static void calcM3i(const double *w, const double *dxv, const double *f, double *out) {
      MomentCalc${ci}x${vi}v${basisShortNm}_M3i_P${pi}(w, dxv, f, out);
    }
    
};

|    end
|  end
|end

]]

-- Instantiate template.
local momentHeaderTemplate = xsys.template ( momentHeaderTemplateString )

-- Concatenate various pieces to generate final header.
momentHeaderOut =
   momentHeaderTemplateTopString
   ..
   momentHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMIN = 1, VMAX = 2,
      PMIN = 1, PMAX = 3,
      basisNm = 'G_TENSOR_C',
      basisShortNm = 'Tensor',
   }
   ..
   momentHeaderTemplate {
      CMIN = 1, CMAX = 1,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 3,
      basisNm = 'G_TENSOR_C',
      basisShortNm = 'Tensor',
   }
   ..
   momentHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMIN = 1, VMAX = 2,
      PMIN = 1, PMAX = 3,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   momentHeaderTemplate {
      CMIN = 1, CMAX = 1,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 3,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   momentHeaderTemplate {
      CMIN = 2, CMAX = 2,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   momentHeaderTemplate {
      CMIN = 3, CMAX = 3,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   momentHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMIN = 1, VMAX = 2,
      PMIN = 1, PMAX = 3,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   momentHeaderTemplate {
      CMIN = 1, CMAX = 1,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 3,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   momentHeaderTemplate {
      CMIN = 2, CMAX = 2,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   momentHeaderTemplate {
      CMIN = 3, CMAX = 3,
      VMIN = 3, VMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   momentHeaderTemplateBottomString

-- Write out string.
fileh = io.open(outDir .. outFileName, 'w')
io.output(fileh)
io.write( momentHeaderOut )
io.close(fileh)

