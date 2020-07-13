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
      const double *ql, const double*  qr,
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
      const double *ql, const double* qr,
      double *outl, double *outr);
  };

  // Template specializations provide actual calls




  template <>
  class MaxwellModDecl<1,1,G_MAX_ORDER_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol1xMaxP1(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf1xMax_X_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class MaxwellModDecl<1,2,G_MAX_ORDER_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol1xMaxP2(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf1xMax_X_P2(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class MaxwellModDecl<2,1,G_MAX_ORDER_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol2xMaxP1(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf2xMax_X_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 1:
              return MaxwellSurf2xMax_Y_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class MaxwellModDecl<2,2,G_MAX_ORDER_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol2xMaxP2(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf2xMax_X_P2(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 1:
              return MaxwellSurf2xMax_Y_P2(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };





  template <>
  class MaxwellModDecl<3,1,G_MAX_ORDER_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol3xMaxP1(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf3xMax_X_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 1:
              return MaxwellSurf3xMax_Y_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 2:
              return MaxwellSurf3xMax_Z_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };





  template <>
  class MaxwellModDecl<1,1,G_SERENDIPITY_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol1xSerP1(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf1xSer_X_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class MaxwellModDecl<1,2,G_SERENDIPITY_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol1xSerP2(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf1xSer_X_P2(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class MaxwellModDecl<2,1,G_SERENDIPITY_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol2xSerP1(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf2xSer_X_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 1:
              return MaxwellSurf2xSer_Y_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


  template <>
  class MaxwellModDecl<2,2,G_SERENDIPITY_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol2xSerP2(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf2xSer_X_P2(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 1:
              return MaxwellSurf2xSer_Y_P2(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };





  template <>
  class MaxwellModDecl<3,1,G_SERENDIPITY_C> {
    public:

      static __host__ __device__ double volumeTerm(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out) {
        return MaxwellVol3xSerP1(meq, w, dx, q, out);
      }

      static __host__ __device__ double surfTerm(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr) {

        switch (dir) {
          case 0:
              return MaxwellSurf3xSer_X_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 1:
              return MaxwellSurf3xSer_Y_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
          case 2:
              return MaxwellSurf3xSer_Z_P1(meq, wl, wr, dxl, dxr, tau, ql, qr, outl, outr);
              break;
        }        
        return 0;
      }
  };


}

