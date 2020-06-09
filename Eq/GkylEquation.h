#ifndef GKYL_EQUATION_H
#define GKYL_EQUATION_H

#include <GkylCudaConfig.h>

extern "C"  {
    
    // Function pointer types for volume and surface terms
    typedef double (*volTermFunc_t)(void *self, 
      double *xc, double *dx, int *idx, double *qIn, double *qRhsOut);
    
    typedef double (*surfTermFunc_t)(void *self, int dir, double *cflL, double *cflR,
      double *xcL, double *xcR, double *dxL, double *dxR,
      double maxsOld, int* idxL, int *idxR,
      double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR);

    typedef double (*boundarySurfTermFunc_t)(void *self, int dir, double *cflL, double *cflR,
      double *xcL, double *xcR, double *dxL, double *dxR,
      double maxsOld, int* idxL, int *idxR,
      double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR);

    typedef struct {
        // Pointer to specific equation object
        void *equation;
        // pointer to equation volTerm function
        volTermFunc_t equationVolTerm;
        // pointer to equation surfTerm function
        surfTermFunc_t equationSurfTerm;
         // pointer to equation boundarySurfTerm function
        boundarySurfTermFunc_t equationBoundarySurfTerm;
        
        __host__ __device__ double volTerm(double *xc, double *dx, int *idx, double *qIn, double *qRhsOut) {
          return equationVolTerm(equation, xc, dx, idx, qIn, qRhsOut);
        }

        __host__ __device__ double surfTerm(int dir, double *cflL, double *cflR,
          double *xcL, double *xcR, double *dxL, double *dxR,
          double maxsOld, int* idxL, int *idxR,
          double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {
          return equationSurfTerm(equation, dir, cflL, cflR,
            xcL, xcR, dxL, dxR,
            maxsOld, idxL, idxR,
            qInL, qInR, qRhsOutL, qRhsOutR);
        }

        __host__ __device__ double boundarySurfTerm(int dir, double *cflL, double *cflR,
          double *xcL, double *xcR, double *dxL, double *dxR,
          double maxsOld, int* idxL, int *idxR,
          double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {
          return equationBoundarySurfTerm(equation, dir, cflL, cflR,
            xcL, xcR, dxL, dxR,
            maxsOld, idxL, idxR,
            qInL, qInR, qRhsOutL, qRhsOutR);
        }
    } GkylEquation_t;
}

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
