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
    }
    ~Vlasov() = default;

    __host__ __device__ int getCdim() { return cdim; }

    __host__ __device__ GkylCartField_t* getEmField() { return emField; }
  
    __host__ __device__ void setAuxFields(GkylCartField_t *emField);
  
    __host__ __device__ double volTerm(const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut);
    __device__ double volTerm_shared(const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut);

    __host__ __device__ double surfTerm(const int dir, 
                    const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
                    const double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
                    const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR);

    __host__ __device__ inline double boundarySurfTerm(const int dir, 
                    const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
                    const double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
                    const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR) {return 0;};

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
    // hard-code a specific dimensionality here for now (polymorphism abstraction TBD)
    // template parameters are 
    // VlasovModDecl<cdim, vdim, polyOrder, basisType>
    // this needs to match the basis used to initialize the equation object in Lua 
    VlasovModDecl<2,3,1,G_SERENDIPITY_C> kernel;
  };
}
#endif
