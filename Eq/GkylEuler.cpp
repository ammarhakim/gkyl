#include <GkylEuler.h>

namespace Gkyl {
  // C wrappers to member functions
  void* new_Euler() {
    Euler *euler = new Euler();
    return reinterpret_cast<void*>(euler);
  }

  void* new_Euler_onDevice(Euler *v) {
    Euler *d_v;
    cudaMalloc((void**)&d_v, sizeof(Euler));
    cudaMemcpy(d_v, v, sizeof(Euler), cudaMemcpyHostToDevice);
    return reinterpret_cast<void*>(d_v);
  }

  __host__ __device__ int Euler::getNumWaves() { return _numWaves; }

  __host__ __device__ int Euler::getNumEquations() { return _numEquations; }

  __host__ __device__ double Euler::gasGamma() { return _gasGamma; }

  __host__ __device__ double Euler::pressure(const double *q) {
    return (_gasGamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
  }
}
