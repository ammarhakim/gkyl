// Gkyl ------------------------------------------------------------------------
//
// Maxwell equation object
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylMaxwell.h>
#include <GkylCudaFuncs.h>
#include <MaxwellModDecl.h>

// cuda includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_types.h>

Maxwell_volumeTerm_t Maxwell_getVolumeTerm(int cdim, int polyOrder, int basisType);
Maxwell_surfTerm_t Maxwell_getSurfTerm(int cdim, int polyOrder, int basisType);

__device__ double Maxwell_volTerm(void *self, 
  double *xc, double *dx, int *idx, double *qIn, double *qRhsOut) {
  
  GkylMaxwell_t *eqn = (GkylMaxwell_t *) self;
  return eqn->volumeTerm(&eqn->mdata, xc, dx, qIn, qRhsOut);
}
__device__ volTermFunc_t p_Maxwell_volTerm = &Maxwell_volTerm;

__device__ double Maxwell_surfTerm(void *self, int dir,
  double *xcL, double *xcR, double *dxL, double *dxR,
  double maxsOld, int* idxL, int *idxR,
  double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {

  GkylMaxwell_t *eqn = (GkylMaxwell_t *) self;
  return eqn->surfTerm(dir, &eqn->mdata, xcL, xcR, dxL, dxR, eqn->tau, qInL, qInR, qRhsOutL, qRhsOutR);
}
__device__ surfTermFunc_t p_Maxwell_surfTerm = &Maxwell_surfTerm;

GkylEquation_t*
new_MaxwellOnDevice(unsigned cdim, unsigned polyOrder, unsigned basisType, MaxwellEq_t *mdata,
  double tau) {
  GkylEquation_t *eqn = new GkylEquation_t;
  GkylMaxwell_t *maxwellEqn = new GkylMaxwell_t;

  maxwellEqn->cdim = cdim;
  maxwellEqn->polyOrder = polyOrder;
  maxwellEqn->basisType = basisType;

  maxwellEqn->mdata.c = mdata->c;
  maxwellEqn->mdata.chi = mdata->chi;
  maxwellEqn->mdata.gamma = mdata->gamma;

  maxwellEqn->tau = tau;

  maxwellEqn->volumeTerm = Maxwell_getVolumeTerm(cdim, polyOrder, basisType);
  maxwellEqn->surfTerm = Maxwell_getSurfTerm(cdim, polyOrder, basisType);

  // setup equation object with Maxwell data and function pointers
  eqn->equation = Gkyl::CudaUtils<GkylMaxwell_t>::allocAndCopyToDevice(maxwellEqn);
  
  // copy functions from device and set them appropriately
  auto err1 = cudaMemcpyFromSymbol(&eqn->equationVolTerm, 
    p_Maxwell_volTerm, sizeof(volTermFunc_t));
  auto err2 = cudaMemcpyFromSymbol(&eqn->equationSurfTerm, 
    p_Maxwell_surfTerm, sizeof(surfTermFunc_t));

  return Gkyl::CudaUtils<GkylEquation_t>::allocAndCopyToDevice(eqn);
}
