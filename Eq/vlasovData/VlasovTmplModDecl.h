#ifndef VLASOV_TMPL_MOD_DECL_H
#define VLASOV_TMPL_MOD_DECL_H

#include <VlasovModDecl.h>
#include <GkBasisTypes.h>

namespace Gkyl {

  // Provides a templated wrapper around the low level C-style kernels
  // so they can be called more systematically from C++ code
  template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
  class VlasovModDecl {
    public:
      /**
       * Volume streaming term
       */
      double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out);

      /**
       * Surface streaming term
       */
      void surfStreamTerm(unsigned dir, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr);
      
  };

  // Template specializations provide actual calls
  template <>
  class VlasovModDecl<1,1,1,G_MAX_ORDER_C> {
    public:
      double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vMaxP1(w, dxv, f, out);
      }
      
      void surfStreamTerm(unsigned dir, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) {
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }
  };
}
#endif // VLASOV_TMPL_MOD_DECL_H
