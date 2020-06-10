local xsys = require "xsys"

local maxwellHeaderTemplateTopString = [[
// Gkyl ------------------------------------------------------------------------
//
// Wrappers for calling into Maxwell kernels
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylCudaConfig.h>
#include <MaxwellModDecl.h>
#include <GkylBasisTypes.h>

extern "C" {

    typedef double (*Maxwell_volumeTerm_t)(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out);
    typedef double (*Maxwell_surfTerm_t)(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      double *ql, const double*  qr,
      double *outl, double *outr);
}

namespace Gkyl {

  // Provides a templated wrapper around the low level C-style kernels
  // so they can be called more systematically from C++ code
  template <unsigned CDIM, unsigned POLYORDER, unsigned BASISTYPE>
  class MaxwellModDecl {
    public:
      /**
       * Volume term
       */
      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out);

      /**
       * Surface terms
       */
      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      double *ql, const double*  qr,
      double *outl, double *outr);
  };

  // Template specializations provide actual calls

]]

local maxwellHeaderTemplateBottomString = [[
}

]]

local maxwellHeaderTemplateString = [[


|for ci = CMIN, CMAX do
|  for pi = PMIN, PMAX do

  template <>
  class MaxwellModDecl<${ci},${pi},${basisNm}> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol${ci}x${basisShortNm}P${pi}(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf${ci}x${basisShortNm}_X_P${pi}(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
|if ci > 1 then
          case 1:
              return MaxwellSurf${ci}x${basisShortNm}_Y_P${pi}(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
|end
|if ci > 2 then
          case 2:
              return MaxwellSurf${ci}x${basisShortNm}_Z_P${pi}(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
|end
        }        
        return 0;
      }
  };

|    end
|end

]]

-- instantiate template
local maxwellHeaderTemplate = xsys.template ( maxwellHeaderTemplateString )

-- concatinate various pieces to generate final header
local maxwellHeaderOut =
   maxwellHeaderTemplateTopString
   ..
   maxwellHeaderTemplate {
      CMIN = 1, CMAX = 2,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   maxwellHeaderTemplate {
      CMIN = 3, CMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   maxwellHeaderTemplate {
      CMIN = 1, CMAX = 2,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   maxwellHeaderTemplate {
      CMIN = 3, CMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   maxwellHeaderTemplateBottomString

-- write out string
io.write( maxwellHeaderOut )
