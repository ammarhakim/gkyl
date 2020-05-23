#pragma once

// Gkyl includes
#include <GkBasisTypes.h>
#include <VlasovTmplModDecl.h>

// std includes
#include <map>
#include <string>
#include <iostream>
#include <memory>

namespace Gkyl {
  class KernelDriver {
    public:
      /**
       * @param cdim Configuration space dimensions
       * @param vdim Velocity space dimensions
       * @param polyOrder Polynomial order
       * @param basisType Type of basis functions (see Lib/GkBasisTypes.h)
       */
      KernelDriver(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType)
        : cdim(cdim), vdim(vdim), polyOrder(polyOrder), basisType(basisType) {
      }
      
      virtual ~KernelDriver() { }

      /**
       * Run kernels in cell-based update mode
       *
       * @param n Number of times to run kernels
       */
      virtual void cellUpdate(unsigned n) = 0;

    protected:
      unsigned cdim, vdim, polyOrder, basisType;
  };

  /** Helper class to allow easier creation of kernels */
  
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

        // populate kernel factory map (Why is this needed? The basis
        // idea here is to key off the kernels bases on names created
        // by makeName. This allows easier construction of needed
        // kernel pointer, which would otherwise require a really
        // nasty if- or switch-statement.)
        std::map<std::string, KernelFactoryBase*> kernelFactory = {
          { makeName(1,1,1, G_MAX_ORDER_C), new KernelFactory<1,1,1, G_MAX_ORDER_C>() },
          { makeName(1,1,2, G_MAX_ORDER_C), new KernelFactory<1,1,2, G_MAX_ORDER_C>() },
          { makeName(1,1,3, G_MAX_ORDER_C), new KernelFactory<1,1,3, G_MAX_ORDER_C>() },
          { makeName(1,2,1, G_MAX_ORDER_C), new KernelFactory<1,2,1, G_MAX_ORDER_C>() },
          { makeName(1,2,2, G_MAX_ORDER_C), new KernelFactory<1,2,2, G_MAX_ORDER_C>() },
          { makeName(1,2,3, G_MAX_ORDER_C), new KernelFactory<1,2,3, G_MAX_ORDER_C>() },
          { makeName(1,3,1, G_MAX_ORDER_C), new KernelFactory<1,3,1, G_MAX_ORDER_C>() },
          { makeName(1,3,2, G_MAX_ORDER_C), new KernelFactory<1,3,2, G_MAX_ORDER_C>() },
          { makeName(1,3,3, G_MAX_ORDER_C), new KernelFactory<1,3,3, G_MAX_ORDER_C>() },
          { makeName(2,2,1, G_MAX_ORDER_C), new KernelFactory<2,2,1, G_MAX_ORDER_C>() },
          { makeName(2,2,2, G_MAX_ORDER_C), new KernelFactory<2,2,2, G_MAX_ORDER_C>() },
          { makeName(2,2,3, G_MAX_ORDER_C), new KernelFactory<2,2,3, G_MAX_ORDER_C>() },
          { makeName(2,3,1, G_MAX_ORDER_C), new KernelFactory<2,3,1, G_MAX_ORDER_C>() },
          { makeName(2,3,2, G_MAX_ORDER_C), new KernelFactory<2,3,2, G_MAX_ORDER_C>() },
          { makeName(2,3,3, G_MAX_ORDER_C), new KernelFactory<2,3,3, G_MAX_ORDER_C>() },
          { makeName(3,3,1, G_MAX_ORDER_C), new KernelFactory<3,3,1, G_MAX_ORDER_C>() },
          { makeName(3,3,2, G_MAX_ORDER_C), new KernelFactory<3,3,2, G_MAX_ORDER_C>() },
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

        // now construct pointer to kernel
        kernel = std::unique_ptr<VlasovModDeclBase>
          (kernelFactory[makeName(cdim,vdim,polyOrder,basisType)]->create()) ;
      }

      ~VlasovKernelDriver() = default;

      /**
       * Run kernels in cell-based update mode
       *
       * @param n Number of times to run kernels
       */
      void cellUpdate(unsigned n) {
      }

    private:
      /** Pointer to kernel */
      std::unique_ptr<VlasovModDeclBase> kernel;

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
