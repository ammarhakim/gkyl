// Gkyl ------------------------------------------------------------------------
//
// Wrappers for calling into Vlasov kernels
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylCudaConfig.h>
#include <GkylBasisTypes.h>
#include <VlasovModDecl.h>

namespace Gkyl {

  // Base class so pointers to children can be stored and used
  class VlasovModDeclBase {
    public:

      virtual ~VlasovModDeclBase() {}

      /**
       * Volume streaming term
       */
      __host__ __device__ virtual double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) = 0;

      /**
       * Surface streaming term
       */
      __host__ __device__ virtual void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr, const double *fl, const double *fr,
        double *outl, double *outr) = 0;

      /**
       * Volume term (total surface + force)
       */
      __host__ __device__ virtual double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) = 0;

      /**
       * Surface terms from EM forces
       */
      __host__ __device__ virtual void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) = 0;
  };

  // Provides a templated wrapper around the low level C-style kernels
  // so they can be called more systematically from C++ code
  template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
  class VlasovModDecl : public VlasovModDeclBase {
    public:
      /**
       * Volume streaming term
       */
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out);

      /**
       * Surface streaming term
       */
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr, const double *fl, const double *fr,
        double *outl, double *outr);

      /**
       * Volume term (total surface + force)
       */
      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out);

      /**
       * Surface terms from EM forces
       */
      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr);
  };

  // Template specializations provide actual calls




  template <>
  class VlasovModDecl<1,1,1,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vMaxP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x1vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,1,2,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vMaxP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x1vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,2,1,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vMaxP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x2vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x2vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,2,2,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vMaxP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x2vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x2vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,3,1,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vMaxP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vMax_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x3vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x3vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag1x3vMax_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,3,2,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vMaxP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vMax_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x3vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x3vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag1x3vMax_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,2,1,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vMaxP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x2vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x2vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,2,2,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vMaxP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x2vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x2vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,3,1,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vMaxP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x3vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x3vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag2x3vMax_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,3,2,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vMaxP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x3vMax_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x3vMax_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag2x3vMax_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };





  template <>
  class VlasovModDecl<3,3,1,G_MAX_ORDER_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream3x3vMaxP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol3x3vMaxP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag3x3vMax_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag3x3vMax_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag3x3vMax_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };





  template <>
  class VlasovModDecl<1,1,1,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vSerP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x1vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,1,2,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x1vSerP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x1vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x1vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x1vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,2,1,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vSerP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x2vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x2vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,2,2,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x2vSerP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x2vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x2vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x2vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x2vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,3,1,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vSerP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vSer_X_P1(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x3vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x3vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag1x3vSer_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<1,3,2,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream1x3vSerP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double *fl, const double *fr,
        double *outl, double *outr) {
        
        switch (dir) {
          case 0:
              VlasovSurfStream1x3vSer_X_P2(wl, wr, dxvl, dxvr, fl, fr, outl, outr);
              break;
        }
      }

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol1x3vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag1x3vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag1x3vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag1x3vSer_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,2,1,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vSerP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x2vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x2vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,2,2,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x2vSerP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x2vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x2vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x2vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,3,1,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vSerP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x3vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x3vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag2x3vSer_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


  template <>
  class VlasovModDecl<2,3,2,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream2x3vSerP2(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol2x3vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag2x3vSer_VX_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag2x3vSer_VY_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag2x3vSer_VZ_P2(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };





  template <>
  class VlasovModDecl<3,3,1,G_SERENDIPITY_C> : public VlasovModDeclBase {
    public:
      __host__ __device__ double volumeStreamTerm(const double *w, const double *dxv, const double *f, double *out) {
        return VlasovVolStream3x3vSerP1(w, dxv, f, out);
      }
      
      __host__ __device__ void surfStreamTerm(unsigned dir, const double *wl, const double *wr,
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

      __host__ __device__ double volumeTerm(const double *w, const double *dxv, const double *E, const double *f, double *out) {
        return VlasovVol3x3vSerP1(w, dxv, E, f, out);
      }

      __host__ __device__ void surfElcMagTerm(unsigned dir, const double *wl, const double *wr,
        const double *dxvl, const double *dxvr,
        const double amax, const double *E, const
        double *fl, const double *fr,
        double *outl, double *outr) {

        switch (dir) {
          case 0:
              VlasovSurfElcMag3x3vSer_VX_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 1:
              VlasovSurfElcMag3x3vSer_VY_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
          case 2:
              VlasovSurfElcMag3x3vSer_VZ_P1(wl, wr, dxvl, dxvr, amax, E, fl, fr, outl, outr);
              break;
        }        
      }
  };


}

