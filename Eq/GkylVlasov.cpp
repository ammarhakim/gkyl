#include <GkylVlasov.h>

__global__ void setAuxFieldsOnDevice(Gkyl::Vlasov *v, GkylCartField_t* emField);

namespace Gkyl {

  __host__ __device__ double Vlasov::volTerm(double *xc, double *dx, int *idx, double *qIn, double *qRhsOut) {
    Gkyl::GenIndexer emIdxr = emField->genIndexer();
    double *em = emField->getDataPtrAt(emIdxr.index(idx));
    double amid = kernel.volumeTerm(xc, dx, em, qIn, qRhsOut);
    return amid;
  }

  __host__ __device__ double Vlasov::surfTerm(int dir, double *cflL, double *cflR,
                  double *xcL, double *xcR, double *dxL, double *dxR,
                  double maxsOld, int *idxL, int *idxR,
                  double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR) {
    double amax = 0.0;
    if(dir < cdim) {
      kernel.surfStreamTerm(dir, xcL, xcR, dxL, dxR, qInL, qInR, qRhsOutL, qRhsOutR);
    } else if(hasForceTerm) {
      Gkyl::GenIndexer emIdxr = emField->genIndexer();
      double *em = emField->getDataPtrAt(emIdxr.index(idxL));
      amax = kernel.surfElcMagTerm(dir-cdim, xcL, xcR, dxL, dxR, maxsOld, em, qInL, qInR, qRhsOutL, qRhsOutR);
    }
    return amax;
  }

  __host__ __device__ void Vlasov::setAuxFields(GkylCartField_t *em) {
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
