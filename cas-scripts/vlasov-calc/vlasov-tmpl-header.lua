local xsys = require "xsys"

vlasovHeaderTemplateString = [[

#ifndef VLASOV_TMPL_MOD_DECL_H
#define VLASOV_TMPL_MOD_DECL_H

#include <VlasovModDecl.h>
#include <GkBasisTypes.h>

namespace Gkyl {

  // Base class so pointers to children can be stored and used
  class VlasovModDeclBase {
    public:
      /**
       * Volume streaming term
       */
      virtual double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out);

      /**
       * Surface streaming term
       */
      virtual void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr, const double *fl, const double *fr,
        double *outl, double *outr);

      /**
       * Volume term (total surface + force)
       */
      double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out);

      /**
       * Surface terms from EM forces
       */
      virtual void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr);
  };

  // Provides a templated wrapper around the low level C-style kernels
  // so they can be called more systematically from C++ code
  template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
  class VlasovModDecl : public VlasovModDeclBase {
    public:
      /**
       * Volume streaming term
       */
      double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out);

      /**
       * Surface streaming term
       */
      void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr, const double *fl, const double *fr,
        double *outl, double *outr);

      /**
       * Volume term (total surface + force)
       */
      double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out);

      /**
       * Surface terms from EM forces
       */
      void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr);
  };

  // Template specializations provide actual calls

|for ci = 1, CDIM do
| for vi = ci, VDIM do
|  for pi = 1, PMAX do
|    for bni = 1, #basisNm do

  template <>
  class VlasovModDecl<${ci},${vi},${pi},${basisNm[bni]}> : public VlasovModDeclBase {
    public:
      double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream${ci}x${vi}v${basisShortNm[bni]}P${pi}(w, dxv, f, out);
      }
      
      void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream${ci}x${vi}v${basisShortNm[bni]}_X_P${pi}(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
|if ci>1 then
          case 1:
              VlasovSurfStream${ci}x${vi}v${basisShortNm[bni]}_Y_P${pi}(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
|end
|if ci>2 then
          case 2:
              VlasovSurfStream${ci}x${vi}v${basisShortNm[bni]}_Z_P${pi}(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
|end
        }
      }

      double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol${ci}x${vi}v${basisShortNm[bni]}P1(w, dxv, E, f, out);
      }

      void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag${ci}x${vi}v${basisShortNm[bni]}_VX_P${pi}(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
|if vi > 1 then
          case 1:
              VlasovSurfElcMag${ci}x${vi}v${basisShortNm[bni]}_VY_P${pi}(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
|end
|if vi > 2 then
          case 2:
              VlasovSurfElcMag${ci}x${vi}v${basisShortNm[bni]}_VZ_P${pi}(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
|end
        }        
      }
  };

|      end
|    end
|  end
|end
}
#endif // VLASOV_TMPL_MOD_DECL_H

]]

-- instantiate template
local vlasovHeaderTemplate = xsys.template ( vlasovHeaderTemplateString )
-- environment to pass to template 
local vlasovEnv = {
   CDIM = 2, VDIM = 3, PMAX = 2,
   basisNm = {'G_SERENDIPITY_C'},
   basisShortNm = {'Ser'},
}

--   basisNm = {'G_MAX_ORDER_C', 'G_SERENDIPITY_C'},
--   basisShortNm = {'Max', 'Ser'},


-- write out string
io.write( vlasovHeaderTemplate (vlasovEnv) )
			   
