#include <GkylVlasov.h>

// cuda includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_types.h>

__global__ void setAuxFieldsOnDevice(Gkyl::Vlasov *v, GkylCartField_t* emField);

namespace Gkyl {

  __host__ __device__ double Vlasov::volTerm(const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut) {
    Gkyl::GenIndexer emIdxr = emField->genIndexer();
    const double *em = emField->getDataPtrAt(emIdxr.index(idx));
    double amid = kernel.volumeTerm(xc, dx, em, qIn, qRhsOut);
    return amid;
  }

  __device__ double Vlasov::volTerm_shared(const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut) {
    Gkyl::GenIndexer emIdxr = emField->genIndexer();
    //double *em = shared;
    //for(int j=0; j<numComponents; j++) {
    //  // linearIdxC can have jumps, just needs to be contiguous for max(32, numComponents) elements
    //  em[j+numComponents*threadIdx.x] = emField->_data[emIdxr.index(idx) + blockIdx.x*j];
    //}
    const double *em = emField->getDataPtrAt(emIdxr.index(idx));
    double amid = kernel.volumeTerm(xc, dx, em, qIn, qRhsOut);
    return amid;
  }

  __host__ __device__ double Vlasov::surfTerm(const int dir, 
                  const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
                  const double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
                  const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR) {
    double amax = 0.0;
    if(dir < cdim) {
      kernel.surfStreamTerm(dir, xcL, xcR, dxL, dxR, qInL, qInR, qRhsOutL, qRhsOutR);
    } else if(hasForceTerm) {
      Gkyl::GenIndexer emIdxr = emField->genIndexer();
      const double *em = emField->getDataPtrAt(emIdxr.index(idxL));
      amax = kernel.surfElcMagTerm(dir-cdim, xcL, xcR, dxL, dxR, maxsOld, em, qInL, qInR, qRhsOutL, qRhsOutR);
    }
    return amax;
  }

  __host__ __device__ void Vlasov::setAuxFields(GkylCartField_t* em) {
    emField = em;
  }

  // C wrappers to member functions
  void* new_Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType, double qbym, bool hasForceTerm) {
    Vlasov *vlasov = new Vlasov(cdim, vdim, polyOrder, basisType, qbym, hasForceTerm);
    return reinterpret_cast<void*>(vlasov);
  }

  void* new_Vlasov_onDevice(Vlasov *v) {
    Vlasov *d_v;
    cudaMalloc((void**)&d_v, sizeof(Vlasov));
    cudaMemcpy(d_v, v, sizeof(Vlasov), cudaMemcpyHostToDevice);
    return reinterpret_cast<void*>(d_v);
  }
  
  void setAuxFields(Vlasov *v, GkylCartField_t *emField) {
    setAuxFieldsOnDevice<<<1, 1>>>(v, emField);
  }
  
  int getCdim(Vlasov *v) {
    return (int) v->getCdim();
  }
}

__global__ void setAuxFieldsOnDevice(Gkyl::Vlasov *v, GkylCartField_t* emField) {
  v->setAuxFields(emField);
}

//--------- struct-base stuff (remove this line and other code when
//--------- this stuff works and is hooked into HyperDistCont


__global__ void VlasovEquation_setAuxFieldsOnDevice(GkylEquation_t *eqn, GkylCartField_t* emField) {
  GkylVlasovEquation_t *self = (GkylVlasovEquation_t *) eqn->equation;
  self->emField = emField;
}

__device__ double VlasovEquation_volTerm(void *eqn, 
  double *xc, double *dx, int *idx, double *qIn, double *qRhsOut) {

  GkylVlasovEquation_t *self = (GkylVlasovEquation_t *) eqn;
  Gkyl::GenIndexer emIdxr = self->emField->genIndexer();
  const double *em = self->emField->getDataPtrAt(emIdxr.index(idx));
  double amid = 0; /*kernel.volumeTerm(xc, dx, em, qIn, qRhsOut);*/
  return amid;
}
__device__ volTermFunc_t p_VlasovEquation_volTerm = &VlasovEquation_volTerm;

__device__ double VlasovEquation_surfTerm(void *eqn, int dir, double *cflL, double *cflR,
  double *xcL, double *xcR, double *dxL, double *dxR,
  double maxsOld, int* idxL, int *idxR,
  double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {

  GkylVlasovEquation_t *self = (GkylVlasovEquation_t *) eqn;
  double amax = 0.0;
  if(dir < self->cdim) {
    /*kernel.surfStreamTerm(dir, xcL, xcR, dxL, dxR, qInL, qInR, qRhsOutL, qRhsOutR);*/
  } else if (self->hasForceTerm) {
    Gkyl::GenIndexer emIdxr = self->emField->genIndexer();
    const double *em = self->emField->getDataPtrAt(emIdxr.index(idxL));
    amax = 0; /*kernel.surfElcMagTerm(dir-cdim, xcL, xcR, dxL, dxR, maxsOld, em, qInL, qInR, qRhsOutL, qRhsOutR);*/
  }
  return amax;
}
__device__ surfTermFunc_t p_VlasovEquation_surfTerm = &VlasovEquation_surfTerm;

__device__ double VlasovEquation_boundarySurfTerm(void *eqn, int dir, double *cflL, double *cflR,
  double *xcL, double *xcR, double *dxL, double *dxR,
  double maxsOld, int* idxL, int *idxR,
  double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {

  GkylVlasovEquation_t *self = (GkylVlasovEquation_t *) eqn;
  return 0;
}
__device__ boundarySurfTermFunc_t p_VlasovEquation_boundarySurfTerm = &VlasovEquation_boundarySurfTerm;

void
VlasovEquation_setAuxFields(GkylEquation_t *eqn, GkylCartField_t* em) {
  VlasovEquation_setAuxFieldsOnDevice<<<1,1>>>(eqn, em);
}

GkylEquation_t*
new_VlasovEquationOnDevice(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType,
  double qbym, bool hasForceTerm) {
  GkylEquation_t *eqn = new GkylEquation_t;
  GkylVlasovEquation_t *vlasovEqn = new GkylVlasovEquation_t;

  vlasovEqn->cdim = cdim;
  vlasovEqn->vdim = vdim;
  vlasovEqn->polyOrder = polyOrder;
  vlasovEqn->basisType = basisType;
  vlasovEqn->qbym = qbym;
  vlasovEqn->hasForceTerm = hasForceTerm;

  // setup equation object with Vlasov data and function pointers
  eqn->equation = vlasovEqn;
  
  // copy functions from device and set them appropriately
  auto err1 = cudaMemcpyFromSymbol(&eqn->equationVolTerm, 
    p_VlasovEquation_volTerm, sizeof(volTermFunc_t));
  auto err2 = cudaMemcpyFromSymbol(&eqn->equationSurfTerm, 
    p_VlasovEquation_surfTerm, sizeof(surfTermFunc_t));
  auto err3 = cudaMemcpyFromSymbol(&eqn->equationBoundarySurfTerm, 
    p_VlasovEquation_boundarySurfTerm, sizeof(boundarySurfTermFunc_t));
  
  return eqn;
}
