#pragma once

// Gkyl includes
#include <GkBasisTypes.h>

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
      virtual void cellBasedUpdate(unsigned n) = 0;

      /**
       * Run kernels in face-based update mode
       *
       * @param n Number of times to run kernels
       */
      virtual void faceBasedUpdate(unsigned n) = 0;

    protected:
      unsigned cdim, vdim, polyOrder, basisType;
  };
}
