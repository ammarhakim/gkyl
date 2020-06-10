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
    typedef void (*Vlasov_surfSreamTerm_t)(unsigned dir, const double*  wl, const double*  wr,
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




  template <>
  class VlasovModDecl<1,1,1,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vMaxP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vMaxP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x1vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,1,2,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vMaxP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vMaxP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x1vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,2,1,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vMaxP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vMaxP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x2vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x2vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,2,2,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vMaxP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vMaxP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x2vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x2vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,3,1,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vMaxP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vMaxP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x3vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x3vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag1x3vMax_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,3,2,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vMaxP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vMaxP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x3vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x3vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag1x3vMax_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,2,1,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vMaxP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x2vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x2vMax_Y_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vMaxP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x2vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x2vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,2,2,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vMaxP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x2vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x2vMax_Y_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vMaxP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x2vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x2vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,3,1,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vMaxP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x3vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x3vMax_Y_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vMaxP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x3vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x3vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag2x3vMax_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,3,2,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vMaxP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x3vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x3vMax_Y_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vMaxP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x3vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x3vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag2x3vMax_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };





  template <>
  class VlasovModDecl<3,3,1,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream3x3vMaxP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream3x3vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream3x3vMax_Y_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfStream3x3vMax_Z_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol3x3vMaxP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag3x3vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag3x3vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag3x3vMax_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<3,3,2,G_MAX_ORDER_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream3x3vMaxP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream3x3vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream3x3vMax_Y_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfStream3x3vMax_Z_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol3x3vMaxP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag3x3vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag3x3vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag3x3vMax_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };





  template <>
  class VlasovModDecl<1,1,1,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vSerP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vSerP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x1vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,1,2,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vSerP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vSerP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x1vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,2,1,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vSerP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vSerP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x2vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x2vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,2,2,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vSerP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vSerP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x2vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x2vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,3,1,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vSerP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vSerP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x3vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x3vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag1x3vSer_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<1,3,2,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vSerP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vSerP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag1x3vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag1x3vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag1x3vSer_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,2,1,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vSerP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x2vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x2vSer_Y_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vSerP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x2vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x2vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,2,2,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vSerP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x2vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x2vSer_Y_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vSerP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x2vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x2vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,3,1,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vSerP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x3vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x3vSer_Y_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vSerP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x3vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x3vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag2x3vSer_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class VlasovModDecl<2,3,2,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vSerP2(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream2x3vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream2x3vSer_Y_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vSerP2(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag2x3vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag2x3vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag2x3vSer_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };





  template <>
  class VlasovModDecl<3,3,1,G_SERENDIPITY_C> {
    public:
      static __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream3x3vSerP1(w, dxv, f, out);
      }
      
      static __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream3x3vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfStream3x3vSer_Y_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfStream3x3vSer_Z_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      static __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol3x3vSerP1(w, dxv, E, f, out);
      }

      static __host__ __device__ double surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              return VlasovSurfElcMag3x3vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              return VlasovSurfElcMag3x3vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              return VlasovSurfElcMag3x3vSer_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
        return 0;
      }
  };


}

