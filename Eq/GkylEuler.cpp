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

  __host__ __device__ void Euler::qFluctuationsLax(
      const unsigned dir, const double *ql, const double *qr,
      const double *waves, const double *s, double *amdq,
      double *apdq) {
    const int *d = dirShuffle[dir];

   // Arrange fluctuations to sum to jump in physical flux.
   auto fl = _fl;
   auto fr = _fr;

   // Left fluxes.
   auto pr = pressure(ql);
   auto u = ql[d[1]]/ql[0];
   fl[0] = ql[d[1]]; // rho*u
   fl[d[1]] = ql[d[1]]*u + pr; // rho*u*u + p
   fl[d[2]] = ql[d[2]]*u; // rho*v*u
   fl[d[3]] = ql[d[3]]*u; // rho*w*u
   fl[4] = (ql[4]+pr)*u; // (E+p)*u
   auto absMaxsl = abs(u)+sqrt(_gasGamma*pr/ql[1]);

   // Right fluxes.
   pr = pressure(qr);
   u = qr[d[1]]/qr[0];
   fr[0] = qr[d[1]]; // rho*u
   fr[d[1]] = qr[d[1]]*u + pr; // rho*u*u + p
   fr[d[2]] = qr[d[2]]*u; // rho*v*u
   fr[d[3]] = qr[d[3]]*u; // rho*w*u
   fr[4] = (qr[4]+pr)*u; // (E+p)*u
   auto absMaxsr =abs(u)+sqrt(_gasGamma*pr/qr[1]);

   auto absMaxs = fmaxf(absMaxsl, absMaxsr);

   // Left going fluctuations.
   amdq[0] = 0.5*(fr[0]-fl[0] - absMaxs*(qr[0]-ql[0]));
   amdq[d[1]] = 0.5*(fr[d[1]]-fl[d[1]] - absMaxs*(qr[d[1]]-ql[d[1]]));
   amdq[d[2]] = 0.5*(fr[d[2]]-fl[d[2]] - absMaxs*(qr[d[2]]-ql[d[2]]));
   amdq[d[3]] = 0.5*(fr[d[3]]-fl[d[3]] - absMaxs*(qr[d[3]]-ql[d[3]]));
   amdq[4] = 0.5*(fr[4]-fl[4] - absMaxs*(qr[4]-ql[4]));

   // Right going fluctuations.
   apdq[0]  = 0.5*(fr[0]-fl[0] + absMaxs*(qr[0]-ql[0]));
   apdq[d[1]] = 0.5*(fr[d[1]]-fl[d[1]] + absMaxs*(qr[d[1]]-ql[d[1]]));
   apdq[d[2]] = 0.5*(fr[d[2]]-fl[d[2]] + absMaxs*(qr[d[2]]-ql[d[2]]));
   apdq[d[3]] = 0.5*(fr[d[3]]-fl[d[3]] + absMaxs*(qr[d[3]]-ql[d[3]]));
   apdq[4] = 0.5*(fr[4]-fl[4] + absMaxs*(qr[4]-ql[4]));
  }
}
