#ifndef GKYL_VLASOV_H
#define GKYL_VLASOV_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylCartField.h>
#include <GkylEquation.h>
#include <VlasovTmplModDecl.h>
#include <GkylBasisTypes.h>

// std includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Gkyl {

  class Vlasov;
  
  /* C wrappers to member functions, so that they can be called from Lua */
  extern "C" {
    void* new_Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType, double qbym, bool hasForceTerm);
    void* new_Vlasov_onDevice(Vlasov *v);
    int getCdim(Vlasov *v);
    void setAuxFields(Vlasov *eq, GkylCartField_t *emField);
  }
  
  class Vlasov {
   public:
    __host__ __device__ Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType, double qbym, bool hasForceTerm) 
     : cdim(cdim), vdim(vdim), polyOrder(polyOrder), basisType(basisType), qbym(qbym), hasForceTerm(hasForceTerm) {
      // populate kernel factory map (Why is this needed? Basic idea
      // here is to key off kernels based on names created by
      // makeName. This allows easier construction of needed kernel
      // pointer, which would otherwise require a really nasty if-
      // or switch-statement.)
//      std::map<int, KernelFactoryBase*> kernelFactory = {
//        { makeName(1,1,1, G_MAX_ORDER_C), new KernelFactory<1,1,1, G_MAX_ORDER_C>() },
//        { makeName(1,1,2, G_MAX_ORDER_C), new KernelFactory<1,1,2, G_MAX_ORDER_C>() },
//        { makeName(1,2,1, G_MAX_ORDER_C), new KernelFactory<1,2,1, G_MAX_ORDER_C>() },
//        { makeName(1,2,2, G_MAX_ORDER_C), new KernelFactory<1,2,2, G_MAX_ORDER_C>() },
//        { makeName(1,3,1, G_MAX_ORDER_C), new KernelFactory<1,3,1, G_MAX_ORDER_C>() },
//        { makeName(1,3,2, G_MAX_ORDER_C), new KernelFactory<1,3,2, G_MAX_ORDER_C>() },
//        { makeName(2,2,1, G_MAX_ORDER_C), new KernelFactory<2,2,1, G_MAX_ORDER_C>() },
//        { makeName(2,2,2, G_MAX_ORDER_C), new KernelFactory<2,2,2, G_MAX_ORDER_C>() },
//        { makeName(2,3,1, G_MAX_ORDER_C), new KernelFactory<2,3,1, G_MAX_ORDER_C>() },
//        { makeName(2,3,2, G_MAX_ORDER_C), new KernelFactory<2,3,2, G_MAX_ORDER_C>() },
//        { makeName(3,3,1, G_MAX_ORDER_C), new KernelFactory<3,3,1, G_MAX_ORDER_C>() },
//        { makeName(1,1,1, G_SERENDIPITY_C), new KernelFactory<1,1,1, G_SERENDIPITY_C>() },
//        { makeName(1,1,2, G_SERENDIPITY_C), new KernelFactory<1,1,2, G_SERENDIPITY_C>() },
//        { makeName(1,2,1, G_SERENDIPITY_C), new KernelFactory<1,2,1, G_SERENDIPITY_C>() },
//        { makeName(1,2,2, G_SERENDIPITY_C), new KernelFactory<1,2,2, G_SERENDIPITY_C>() },
//        { makeName(1,3,1, G_SERENDIPITY_C), new KernelFactory<1,3,1, G_SERENDIPITY_C>() },
//        { makeName(1,3,2, G_SERENDIPITY_C), new KernelFactory<1,3,2, G_SERENDIPITY_C>() },
//        { makeName(2,2,1, G_SERENDIPITY_C), new KernelFactory<2,2,1, G_SERENDIPITY_C>() },
//        { makeName(2,2,2, G_SERENDIPITY_C), new KernelFactory<2,2,2, G_SERENDIPITY_C>() },
//        { makeName(2,3,1, G_SERENDIPITY_C), new KernelFactory<2,3,1, G_SERENDIPITY_C>() },
//        { makeName(2,3,2, G_SERENDIPITY_C), new KernelFactory<2,3,2, G_SERENDIPITY_C>() },
//        { makeName(3,3,1, G_SERENDIPITY_C), new KernelFactory<3,3,1, G_SERENDIPITY_C>() },          
//      };
  
      // now construct pointer to kernel (above horrid mess is to
      // enable the following one-line construction of the kernel
      // pointer)
//      kernel = (VlasovModDeclBase*)
//        (kernelFactory[makeName(cdim,vdim,polyOrder,basisType)]->create()) ;
  
    }
    ~Vlasov() = default;

    __host__ __device__ int getCdim() { return cdim; }

    __host__ __device__ GkylCartField_t* getEmField() { return emField; }
  
    __host__ __device__ void setAuxFields(GkylCartField_t *emField);
  
    __host__ __device__ double volTerm(double *xc, double *dx, int *idx, double *qIn, double *qRhsOut);

    __host__ __device__  double surfTerm(int dir, double *cflL, double *cflR,
                    double *xcL, double *xcR, double *dxL, double *dxR,
                    double maxsOld, int *idxL, int *idxR,
                    double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR);

    __host__ __device__ inline double boundarySurfTerm(int dir, double *cflL, double *cflR,
                            double *xcL, double *xcR, double *dxL, double *dxR,
                            double maxsOld, int *idxL, int *idxR,
                            double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {return 0;};
  
   private:
    /* dimension and basis parameters */
    const unsigned cdim;
    const unsigned vdim;
    const unsigned polyOrder;
    const unsigned basisType;
  
    /* species parameters */
    const double qbym;
    const bool hasForceTerm;
   
    /* pointers to fields */
    GkylCartField_t *emField;
  
    /** Pointer to kernel */
    VlasovModDecl<1,2,2,G_SERENDIPITY_C> kernel;
  
    ///** Private classes to allow making kernel objects */
    //class KernelFactoryBase {
    //  public:
    //    virtual __host__ __device__ Gkyl::VlasovModDeclBase* create() = 0;
    //};
    //
    //template <unsigned CDIM, unsigned VDIM, unsigned POLYORDER, unsigned BASISTYPE>
    //class KernelFactory : public KernelFactoryBase {
    //  public:
    //    __host__ __device__ Gkyl::VlasovModDeclBase* create() {
    //      return new Gkyl::VlasovModDecl<CDIM, VDIM, POLYORDER, BASISTYPE>();
    //    }
    //};
  
    ///** Make a unique name */
    //__host__ __device__ int makeName(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType) {
    //  
    //  return cdim+10*vdim+100*polyOrder+1000*basisType;
    //}
  };
}
#endif
