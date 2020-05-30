#ifndef GKYL_EQUATION_H
#define GKYL_EQUATION_H

#include <GkylCudaConfig.h>

class Equation {
 public:
  Equation() {};
  virtual ~Equation() {};

  virtual __host__ __device__ double volTerm(double *xc, double *dx, int *idx, double *qIn, double *qRhsOut)=0;
  virtual __host__ __device__ double surfTerm(int dir, double *cflL, double *cflR,
                           double *xcL, double *xcR, double *dxL, double *dxR,
                           double maxsOld, int* idxL, int *idxR,
                           double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR)=0;
  virtual __host__ __device__ double boundarySurfTerm(int dir, double *cflL, double *cflR,
                          double *xcL, double *xcR, double *dxL, double *dxR,
                          double maxsOld, int* idxL, int *idxR,
                          double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR)=0;

 private:

};

#endif // GKYL_EQUATION_H
