// Gkyl ------------------------------------------------------------------------
//
// Vlasov equation object
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylVlasov.h>
#include <GkylCudaFuncs.h>
#include <VlasovModDecl.h>

// cuda includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_types.h>

Vlasov_volumeStreamTerm_t Vlasov_getVolumeStreamTerm(int cdim, int vdim, int polyOrder, int basisType);
Vlasov_surfStreamTerm_t Vlasov_getSurfStreamTerm(int cdim, int vdim, int polyOrder, int basisType);
Vlasov_volumeTerm_t Vlasov_getVolumeTerm(int cdim, int vdim, int polyOrder, int basisType);
Vlasov_surfElcMagTerm_t Vlasov_getSurfElcMagTerm(int cdim, int vdim, int polyOrder, int basisType);

__global__ void Vlasov_setAuxFieldsOnDevice(GkylEquation_t *eqn, GkylCartField_t* emField) {
  GkylVlasov_t *self = (GkylVlasov_t *) eqn->equation;
  self->emField = emField;
}

void
Vlasov_setAuxFields(GkylEquation_t *eqn, GkylCartField_t* em) {
  Vlasov_setAuxFieldsOnDevice<<<1,1>>>(eqn, em);
}

__device__ double Vlasov_volTerm(const void* __restrict__ eqn, 
  const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut) {

  const GkylVlasov_t *self = (GkylVlasov_t *) eqn;
  Gkyl::GenIndexer emIdxr = self->emField->genIndexer();
  const double *em = self->emField->getDataPtrAt(emIdxr.index(idx));
  double cflFreq = self->volumeTerm(xc, dx, em, qIn, qRhsOut);
  return cflFreq;
}
__device__ volTermFunc_t p_Vlasov_volTerm = &Vlasov_volTerm;

__device__ double Vlasov_surfTerm(const void* __restrict__ eqn, int dir,
  const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
  const double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
  const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR) {

  const GkylVlasov_t *self = (GkylVlasov_t *) eqn;
  double amax = 0.0;
  if (dir < self->cdim) {
    self->surfStreamTerm(dir, xcL, xcR, dxL, dxR, qInL, qInR, qRhsOutL, qRhsOutR);
  } else if (self->hasForceTerm) {
    Gkyl::GenIndexer emIdxr = self->emField->genIndexer();
    const double *em = self->emField->getDataPtrAt(emIdxr.index(idxL));
    amax = self->surfElcMagTerm(dir-self->cdim, xcL, xcR, dxL, dxR, maxsOld, em, qInL, qInR, qRhsOutL, qRhsOutR);
  }
  return amax;
}
__device__ surfTermFunc_t p_Vlasov_surfTerm = &Vlasov_surfTerm;

__device__ double Vlasov_boundarySurfTerm(const void* __restrict__ eqn, int dir,
  const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
  const double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
  const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR) {

  return 0;
}
__device__ boundarySurfTermFunc_t p_Vlasov_boundarySurfTerm = &Vlasov_boundarySurfTerm;

GkylEquation_t*
new_VlasovOnDevice(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType,
  double qbym, bool hasForceTerm) {
  GkylEquation_t *eqn = new GkylEquation_t;
  GkylVlasov_t *vlasovEqn = new GkylVlasov_t;

  vlasovEqn->cdim = cdim;
  vlasovEqn->vdim = vdim;
  vlasovEqn->polyOrder = polyOrder;
  vlasovEqn->basisType = basisType;
  vlasovEqn->qbym = qbym;
  vlasovEqn->hasForceTerm = hasForceTerm;

  vlasovEqn->volumeStreamTerm = Vlasov_getVolumeStreamTerm(cdim, vdim, polyOrder, basisType);
  vlasovEqn->surfStreamTerm = Vlasov_getSurfStreamTerm(cdim, vdim, polyOrder, basisType);
  vlasovEqn->volumeTerm = Vlasov_getVolumeTerm(cdim, vdim, polyOrder, basisType);
  vlasovEqn->surfElcMagTerm = Vlasov_getSurfElcMagTerm(cdim, vdim, polyOrder, basisType);

  // setup equation object with Vlasov data and function pointers
  eqn->equation = Gkyl::CudaUtils<GkylVlasov_t>::allocAndCopyToDevice(vlasovEqn);
  
  // copy functions from device and set them appropriately
  auto err1 = cudaMemcpyFromSymbol(&eqn->equationVolTerm, 
    p_Vlasov_volTerm, sizeof(volTermFunc_t));
  auto err2 = cudaMemcpyFromSymbol(&eqn->equationSurfTerm, 
    p_Vlasov_surfTerm, sizeof(surfTermFunc_t));
  auto err3 = cudaMemcpyFromSymbol(&eqn->equationBoundarySurfTerm, 
    p_Vlasov_boundarySurfTerm, sizeof(boundarySurfTermFunc_t));

  return Gkyl::CudaUtils<GkylEquation_t>::allocAndCopyToDevice(eqn);
}
