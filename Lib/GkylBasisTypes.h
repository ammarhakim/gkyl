// Gkyl ------------------------------------------------------------------------
//
// Basis type declarations
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

// std includes
#include <string>

namespace Gkyl {

  // Basis function types
  enum BasisTypes {
    G_MAX_ORDER_C = 1, // Max-order
    G_SERENDIPITY_C = 2, // serendipity 
    G_TENSOR_C = 3, // tensor-product
  };

  class BasisInfo {
    public:
      /**
       * Number of basis in this set.
       */
      virtual unsigned numBasis() const = 0;

      /**
       * Basis function identifier string
       */
      virtual std::string id() const = 0;
  };

  class CartModalMaxOrderBasisInfo : public BasisInfo {
    public:
      CartModalMaxOrderBasisInfo(unsigned ndim, unsigned polyOrder)
        : ndim(ndim), polyOrder(polyOrder) {
        // number of basis is = (p+d)! / p! d!
        if (ndim == 0) 
          nbasis = 1;
        else if (ndim == 1)
          nbasis = polyOrder+1;
        else if (ndim == 2) 
          nbasis = (polyOrder+2)*(polyOrder+1)/2;
        else if (ndim == 3) 
          nbasis = (polyOrder+3)*(polyOrder+2)*(polyOrder+1)/6;
        else if (ndim == 4) 
          nbasis = (polyOrder+4)*(polyOrder+3)*(polyOrder+2)*(polyOrder+1)/24;
        else if (ndim == 5) 
          nbasis = (polyOrder+5)*(polyOrder+4)*(polyOrder+3)*(polyOrder+2)*(polyOrder+1)/120;
        else if (ndim == 6) 
          nbasis = (polyOrder+6)*(polyOrder+5)*(polyOrder+4)*(polyOrder+3)*(polyOrder+2)*(polyOrder+1)/720;
      }

      /**
       * Number of basis in this set.
       */
      unsigned numBasis() const { return nbasis; }

      /**
       * Basis function identifier string
       */
      std::string id() const { return "maximal-order"; }

    private:
      /** Dimension */
      unsigned ndim;
      /** Poly order */
      unsigned polyOrder;
      /** Number of basis functions */
      unsigned nbasis;
  };

  class CartModalSerendipityBasisInfo : public BasisInfo {
    public:
      CartModalSerendipityBasisInfo(unsigned ndim, unsigned polyOrder)
        : ndim(ndim), polyOrder(polyOrder), nbasis(1) {
        // See Arnold, Awanou, F. Comp Math, 2011
        unsigned numBasis_d2[] = {4, 8, 12};
        unsigned numBasis_d3[] = {8, 20, 32};
        unsigned numBasis_d4[] = {16, 48, 80};
        unsigned numBasis_d5[] = {32, 112, 192};
        unsigned numBasis_d6[] = {64, 256};

        if (polyOrder > 0) {
          if (ndim == 1)
            nbasis = polyOrder+1;
          else if (ndim == 2)
            nbasis = numBasis_d2[polyOrder-1];
          else if (ndim == 3)
            nbasis = numBasis_d3[polyOrder-1];
          else if (ndim == 4)
            nbasis = numBasis_d4[polyOrder-1];
          else if (ndim == 5)
            nbasis = numBasis_d5[polyOrder-1];
          else if (ndim == 6)
            nbasis = numBasis_d6[polyOrder-1];
        }
      }

      /**
       * Number of basis in this set.
       */
      unsigned numBasis() const { return nbasis; }

      /**
       * Basis function identifier string
       */
      std::string id() const { return "serendipity"; }

    private:
      /** Dimension */
      unsigned ndim;
      /** Poly order */
      unsigned polyOrder;
      /** Number of basis functions */
      unsigned nbasis;
  };  
}
