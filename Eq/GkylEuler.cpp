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

  __host__ __device__ void Euler::rpLax(
      const unsigned dir, const double *delta, const double *ql,
      const double *qr, double **waves, double *s) {
    const int *d = dirShuffle[dir];

    double *wv = waves[0];
    wv[0] = qr[0]-ql[0];
    wv[1] = qr[1]-ql[1];
    wv[2] = qr[2]-ql[2];
    wv[3] = qr[3]-ql[3];
    wv[4] = qr[4]-ql[4];

    double pr = pressure(ql);
    double u = ql[d[1]]/ql[0];
    double cs = std::sqrt(_gasGamma*pr/ql[0]);
    double sl = u+cs;

    pr = pressure(qr);
    u = qr[d[1]]/qr[0];
    cs = std::sqrt(_gasGamma*pr/qr[0]);
    double sr = u+cs   ;

    s[1] = 0.5*(sl+sr);
  }
}
