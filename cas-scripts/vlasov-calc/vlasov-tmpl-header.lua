local xsys = require "xsys"

local vlasovHeaderTemplateTopString = [[
// Gkyl ------------------------------------------------------------------------
//
// Wrappers for calling into Vlasov kernels
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylCudaConfig.h>
#include <VlasovModDecl.h>
#include <GkylBasisTypes.h>

extern "C" {
    typedef double (*Vlasov_volumeStreamTerm_t)(const double*  w, const double*  dxv, const double*  f, double *out);
    typedef void (*Vlasov_surfStreamTerm_t)(unsigned dir, const double*  wl, const double*  wr,
      const double*  dxvl, const double*  dxvr, const double*  fl, const double*  fr,
      double *outl, double *outr);

    typedef double (*Vlasov_volumeTerm_t)(const double*  w, const double*  dxv,
      const double*  E, const double*  f, double *out);
    typedef double (*Vlasov_surfElcMagTerm_t)(unsigned dir, const double*  wl, const double*  wr,
      const double*  dxvl, const double*  dxvr,
      const double amax, const double*  E, const
      double *fl, const double*  fr,
      double *outl, double *outr);    
}

namespace Gkyl {

  // Provides a templated wrapper around the low level C-style kernels
  // so they can be called more systematically from C++ code
  template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
  class VlasovModDecl {
    public:
      /**
       * Volume streaming term
       */
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out);

      /**
       * Surface streaming term
       */
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr, const double *fl, const double *fr,
        double *outl, double *outr);

      /**
       * Volume term (total surface + force)
       */
      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out);

      /**
       * Surface terms from EM forces
       */
      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr);
  };

  // Template specializations provide actual calls

]]

local vlasovHeaderTemplateBottomString = [[
}

]]

local vlasovHeaderTemplateString = [[


|for ci = CMIN, CMAX do
| for vi = ci, VMAX do
|  for pi = PMIN, PMAX do

  template <>
  class VlasovModDecl<${ci},${vi},${pi},${basisNm}> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream${ci}x${vi}v${basisShortNm}P${pi}(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream${ci}x${vi}v${basisShortNm}_X_P${pi}(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
|if ci>1 then
          case 1:
              VlasovSurfStream${ci}x${vi}v${basisShortNm}_Y_P${pi}(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
|end
|if ci>2 then
          case 2:
              VlasovSurfStream${ci}x${vi}v${basisShortNm}_Z_P${pi}(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
|end
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol${ci}x${vi}v${basisShortNm}P${pi}(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag${ci}x${vi}v${basisShortNm}_VX_P${pi}(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
|if vi > 1 then
          case 1:
              return VlasovSurfElcMag${ci}x${vi}v${basisShortNm}_VY_P${pi}(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
|end
|if vi > 2 then
          case 2:
              return VlasovSurfElcMag${ci}x${vi}v${basisShortNm}_VZ_P${pi}(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
|end
        }        
        return 0;
      }
  };

|    end
|  end
|end

]]

-- instantiate template
local vlasovHeaderTemplate = xsys.template ( vlasovHeaderTemplateString )

-- concatinate various pieces to generate final header
local vlasovHeaderOut =
   vlasovHeaderTemplateTopString
   ..
   vlasovHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   vlasovHeaderTemplate {
      CMIN = 3, CMAX = 3,
      VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   vlasovHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   vlasovHeaderTemplate {
      CMIN = 3, CMAX = 3,
      VMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   vlasovHeaderTemplateBottomString

-- write out string
io.write( vlasovHeaderOut )
