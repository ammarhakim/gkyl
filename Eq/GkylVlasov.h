#ifndef GKYL_VLASOV_H
#define GKYL_VLASOV_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylRange.h>
#include <GkylRectCart.h>
#include <GkylCartField.h>
#include <GkylEquation.h>

class Vlasov;

/* C wrappers to member functions, so that they can be called from Lua */
extern "C" {
  void* new_Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, const char* basisType, double qbym, bool hasForceTerm);
  void setAuxFields(Vlasov *eq, GkylCartField_t *emField);
}

class Vlasov: public Equation {
 public:
  Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, const char* basisType, double qbym, bool hasForceTerm);
  ~Vlasov();

  void setAuxFields(GkylCartField_t *emField);

  __host__ __device__ double volTerm(double *xc, double *dx, int *idx, double *qIn, double *qRhsOut) {
  }
  __host__ __device__ double surfTerm(int dir, double *cflL, double *cflR,
                  double *xcL, double *xcR, double *dxL, double *dxR,
                  double maxsOld, int* idxL, int *idxR,
                  double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {};
  __host__ __device__ double boundarySurfTerm(int dir, double *cflL, double *cflR,
                          double *xcL, double *xcR, double *dxL, double *dxR,
                          double maxsOld, int* idxL, int *idxR,
                          double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {};

 private:
  /* dimension and basis parameters */
  const unsigned cdim;
  const unsigned vdim;
  const unsigned polyOrder;
  const char* basisType;

  /* species parameters */
  const double qbym;
  const bool hasForceTerm;
  long a;
 
  /* pointers to fields */
  GkylCartField_t *emField;
  GkylCartField_t *dummy;
};

#endif
