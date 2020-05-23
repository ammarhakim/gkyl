#pragma once

// Gkyl includes
#include <GkBasisTypes.h>
#include <GkKernelDriver.h>
#include <VlasovTmplModDecl.h>

// std includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Gkyl {
  
  class VlasovKernelDriver : public KernelDriver {
    public:
      /**
       * @param cdim Configuration space dimensions
       * @param vdim Velocity space dimensions
       * @param polyOrder Polynomial order
       * @param basisType Type of basis functions (see Lib/GkBasisTypes.h)
       */
      VlasovKernelDriver(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType)
        : KernelDriver(cdim, vdim, polyOrder, basisType) {

        // populate kernel factory map (Why is this needed? Basic idea
        // here is to key off kernels based on names created by
        // makeName. This allows easier construction of needed kernel
        // pointer, which would otherwise require a really nasty if-
        // or switch-statement.)
        std::map<std::string, KernelFactoryBase*> kernelFactory = {
          { makeName(1,1,1, G_MAX_ORDER_C), new KernelFactory<1,1,1, G_MAX_ORDER_C>() },
          { makeName(1,1,2, G_MAX_ORDER_C), new KernelFactory<1,1,2, G_MAX_ORDER_C>() },
          { makeName(1,2,1, G_MAX_ORDER_C), new KernelFactory<1,2,1, G_MAX_ORDER_C>() },
          { makeName(1,2,2, G_MAX_ORDER_C), new KernelFactory<1,2,2, G_MAX_ORDER_C>() },
          { makeName(1,3,1, G_MAX_ORDER_C), new KernelFactory<1,3,1, G_MAX_ORDER_C>() },
          { makeName(1,3,2, G_MAX_ORDER_C), new KernelFactory<1,3,2, G_MAX_ORDER_C>() },
          { makeName(2,2,1, G_MAX_ORDER_C), new KernelFactory<2,2,1, G_MAX_ORDER_C>() },
          { makeName(2,2,2, G_MAX_ORDER_C), new KernelFactory<2,2,2, G_MAX_ORDER_C>() },
          { makeName(2,3,1, G_MAX_ORDER_C), new KernelFactory<2,3,1, G_MAX_ORDER_C>() },
          { makeName(2,3,2, G_MAX_ORDER_C), new KernelFactory<2,3,2, G_MAX_ORDER_C>() },
          { makeName(3,3,1, G_MAX_ORDER_C), new KernelFactory<3,3,1, G_MAX_ORDER_C>() },
          { makeName(1,1,1, G_SERENDIPITY_C), new KernelFactory<1,1,1, G_SERENDIPITY_C>() },
          { makeName(1,1,2, G_SERENDIPITY_C), new KernelFactory<1,1,2, G_SERENDIPITY_C>() },
          { makeName(1,2,1, G_SERENDIPITY_C), new KernelFactory<1,2,1, G_SERENDIPITY_C>() },
          { makeName(1,2,2, G_SERENDIPITY_C), new KernelFactory<1,2,2, G_SERENDIPITY_C>() },
          { makeName(1,3,1, G_SERENDIPITY_C), new KernelFactory<1,3,1, G_SERENDIPITY_C>() },
          { makeName(1,3,2, G_SERENDIPITY_C), new KernelFactory<1,3,2, G_SERENDIPITY_C>() },
          { makeName(2,2,1, G_SERENDIPITY_C), new KernelFactory<2,2,1, G_SERENDIPITY_C>() },
          { makeName(2,2,2, G_SERENDIPITY_C), new KernelFactory<2,2,2, G_SERENDIPITY_C>() },
          { makeName(2,3,1, G_SERENDIPITY_C), new KernelFactory<2,3,1, G_SERENDIPITY_C>() },
          { makeName(2,3,2, G_SERENDIPITY_C), new KernelFactory<2,3,2, G_SERENDIPITY_C>() },
          { makeName(3,3,1, G_SERENDIPITY_C), new KernelFactory<3,3,1, G_SERENDIPITY_C>() },          
        };

        // now construct pointer to kernel (above horrid mess is to
        // enable the following one-line construction of the kernel
        // pointer)
        kernel = std::unique_ptr<VlasovModDeclBase>
          (kernelFactory[makeName(cdim,vdim,polyOrder,basisType)]->create()) ;

        // for use in tests, store number of basis functions
        if (basisType == Gkyl::G_MAX_ORDER_C) {
          CartModalMaxOrderBasisInfo confInfo(cdim, polyOrder);
          numConfBasis = confInfo.numBasis();
          CartModalMaxOrderBasisInfo phaseInfo(cdim+vdim, polyOrder);
          numPhaseBasis = phaseInfo.numBasis();          
        }
        else if (basisType == Gkyl::G_SERENDIPITY_C) {
          CartModalSerendipityBasisInfo confInfo(cdim, polyOrder);
          numConfBasis = confInfo.numBasis();
          CartModalSerendipityBasisInfo phaseInfo(cdim+vdim, polyOrder);
          numPhaseBasis = phaseInfo.numBasis();
        }
      }

      ~VlasovKernelDriver() = default;

      /**
       * Run kernels in cell-based update mode
       *
       * @param n Number of times to run kernels
       */
      void cellBasedUpdate(unsigned n) {
        unsigned ndim = cdim+vdim;
        // allocate arrays to pass to kernels
        std::vector<double> w(ndim), dxv(ndim), E(8*numConfBasis), f(numPhaseBasis), out(numPhaseBasis);
        std::vector<double> wl(ndim), wr(ndim), dxvl(ndim), dxvr(ndim), fl(numPhaseBasis), fr(numPhaseBasis);
        std::vector<double> outl(numPhaseBasis), outr(numPhaseBasis);

        for (auto i=0; i<n; ++i)
        {
          // volume term
          kernel->volumeTerm(&w[0], &dxv[0], &E[0], &f[0], &out[0]);
          // surface streaming terms (run surface kernels twice as there are two faces in each direction)
          for (auto d=0; d<cdim; ++d) {
            kernel->surfStreamTerm(d, &wl[0], &wr[0], &dxvl[0], &dxvr[0], &fl[0], &fr[0], &outl[0], &outr[0]);
            kernel->surfStreamTerm(d, &wl[0], &wr[0], &dxvl[0], &dxvr[0], &fl[0], &fr[0], &outl[0], &outr[0]);
          }            
          // force streaming terms (run surface kernels twice as there are two faces in each direction)
          for (auto d=0; d<vdim; ++d) {
            kernel->surfElcMagTerm(d, &wl[0], &wr[0], &dxvl[0], &dxvr[0], 1.0, &E[0],
              &fl[0], &fr[0], &outl[0], &outr[0]);
            kernel->surfElcMagTerm(d, &wl[0], &wr[0], &dxvl[0], &dxvr[0], 1.0, &E[0],
              &fl[0], &fr[0], &outl[0], &outr[0]);
          }
        }
      }

      /**
       * Run kernels in face-based update mode
       *
       * @param n Number of times to run kernels
       */
      void faceBasedUpdate(unsigned n) {
        unsigned ndim = cdim+vdim;
        // allocate arrays to pass to kernels
        std::vector<double> w(ndim), dxv(ndim), E(8*numConfBasis), f(numPhaseBasis), out(numPhaseBasis);
        std::vector<double> wl(ndim), wr(ndim), dxvl(ndim), dxvr(ndim), fl(numPhaseBasis), fr(numPhaseBasis);
        std::vector<double> outl(numPhaseBasis), outr(numPhaseBasis);

        for (auto i=0; i<n; ++i)
        {
          // volume term
          kernel->volumeTerm(&w[0], &dxv[0], &E[0], &f[0], &out[0]);
          // surface streaming terms
          for (auto d=0; d<cdim; ++d)
            kernel->surfStreamTerm(d, &wl[0], &wr[0], &dxvl[0], &dxvr[0], &fl[0], &fr[0], &outl[0], &outr[0]);
          // force streaming terms
          for (auto d=0; d<vdim; ++d)
            kernel->surfElcMagTerm(d, &wl[0], &wr[0], &dxvl[0], &dxvr[0], 1.0, &E[0],
              &fl[0], &fr[0], &outl[0], &outr[0]);
        }
      }      

    private:
      /** Pointer to kernel */
      std::unique_ptr<VlasovModDeclBase> kernel;
      /** Number of conf space basis functions */
      unsigned numConfBasis;
      /** Number of phase-space basis functions */
      unsigned numPhaseBasis;

      /** Private classes to allow making kernel objects */
      class KernelFactoryBase {
        public:
          virtual Gkyl::VlasovModDeclBase* create() = 0;
      };
      
      template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
      class KernelFactory : public KernelFactoryBase {
        public:
          Gkyl::VlasovModDeclBase* create() {
            return new Gkyl::VlasovModDecl<CDIM, VDIM, POLYORDER, BASISTYPE>();
          }
      };

      /** Make a unique name */
      std::string makeName(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType) {
        std::string bc;
        if (basisType == Gkyl::G_MAX_ORDER_C)
          bc = "Max";
        else if (basisType == Gkyl::G_SERENDIPITY_C)
          bc = "Ser";
        else if (basisType == Gkyl::G_TENSOR_C)
          bc = "Ten";
        
        return "b" + std::to_string(cdim) + "x" + std::to_string(vdim) + "v"
          + std::to_string(polyOrder) + "p" + bc;
      }
  };
}
