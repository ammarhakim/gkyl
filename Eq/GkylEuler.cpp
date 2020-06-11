#include <GkylEuler.h>

namespace Gkyl {
  // C wrappers to member functions
  void* new_Euler(const double gasGamma) {
    Euler *euler = new Euler(gasGamma);
    return reinterpret_cast<void*>(euler);
  }

  void* new_Euler_onDevice(Euler *eq) {
    Euler *d_eq;
    cudaMalloc((void**)&d_eq, sizeof(Euler));
    cudaMemcpy(d_eq, eq, sizeof(Euler), cudaMemcpyHostToDevice);
    return reinterpret_cast<void*>(d_eq);
  }

  int numEquations_Euler(Euler *eq) {
    return (int) eq->numEquations();
  }

  void rp_Euler(
      Euler *eq, const int dir, const double *delta, const double *ql,
      const double *qr, double *waves, double *s) {
    eq->rp(dir, delta, ql, qr, waves, s);
  }

  void qFluctuations_Euler(
      Euler *eq, const int dir, const double *ql, const double *qr,
      const double *waves, const double *s, double *amdq, double *apdq) {
    eq->qFluctuations(dir, ql, qr, waves, s, amdq, apdq);
  }

  void flux_Euler(Euler *eq, const int dir, const double *qIn, double *fOut) {
    eq->flux(dir, qIn, fOut);
  }
  
  __host__ __device__ int Euler::numWaves() { return _numWaves; }

  __host__ __device__ int Euler::numEquations() { return _numEquations; }

  __host__ __device__ double Euler::gasGamma() { return _gasGamma; }

  __host__ __device__ double Euler::pressure(const double *q) {
    return (_gasGamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
  }

  __host__ __device__ void Euler::rpLax(
      const int dir, const double *delta, const double *ql,
      const double *qr, double *waves, double *s) {
    const int *d = dirShuffle[dir];

    double *wv = &waves[0];
    wv[0] = qr[0]-ql[0];
    wv[1] = qr[1]-ql[1];
    wv[2] = qr[2]-ql[2];
    wv[3] = qr[3]-ql[3];
    wv[4] = qr[4]-ql[4];

    double pr = pressure(ql);
    double u = ql[d[1]]/ql[0];
    double cs = sqrt(_gasGamma*pr/ql[0]);
    double sl = u+cs;

    pr = pressure(qr);
    u = qr[d[1]]/qr[0];
    cs = sqrt(_gasGamma*pr/qr[0]);
    double sr = u+cs   ;

    s[0] = 0.5*(sl+sr);
  }

  __host__ __device__ void Euler::qFluctuationsLax(
      const int dir, const double *ql, const double *qr,
      const double *waves, const double *s, double *amdq,
      double *apdq) {
    const int *d = dirShuffle[dir];

    // Arrange fluctuations to sum to jump in physical flux.
    double fl[10];
    double fr[10];

    // Left fluxes.
    double pr = pressure(ql);
    double u = ql[d[1]]/ql[0];
    fl[0] = ql[d[1]]; // rho*u
    fl[d[1]] = ql[d[1]]*u + pr; // rho*u*u + p
    fl[d[2]] = ql[d[2]]*u; // rho*v*u
    fl[d[3]] = ql[d[3]]*u; // rho*w*u
    fl[4] = (ql[4]+pr)*u; // (E+p)*u
    double absMaxsl = abs(u)+sqrt(_gasGamma*pr/ql[0]);

    // Right fluxes.
    pr = pressure(qr);
    u = qr[d[1]]/qr[0];
    fr[0] = qr[d[1]]; // rho*u
    fr[d[1]] = qr[d[1]]*u + pr; // rho*u*u + p
    fr[d[2]] = qr[d[2]]*u; // rho*v*u
    fr[d[3]] = qr[d[3]]*u; // rho*w*u
    fr[4] = (qr[4]+pr)*u; // (E+p)*u
    double absMaxsr = abs(u)+sqrt(_gasGamma*pr/qr[0]);

    double absMaxs = max(absMaxsl, absMaxsr);

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

  __host__ __device__ void Euler::rp(
      const int dir, const double *delta, const double *ql,
      const double *qr, double *waves, double *s) {
    rpLax(dir, delta, ql, qr, waves, s);
  }

  __host__ __device__ void Euler::qFluctuations(
      const int dir, const double *ql, const double *qr,
      const double *waves, const double *s, double *amdq,
      double *apdq) {
    qFluctuationsLax(dir, ql, qr, waves, s, amdq, apdq);
  }

  __host__ __device__ void Euler::flux(
      const int dir, const double *qIn, double *fOut) {
    const int *d = dirShuffle[dir];
    double pr = pressure(qIn);
    double u = qIn[d[1]]/qIn[0];
    fOut[0] = qIn[d[1]]; // rho*u.
    fOut[d[1]] = qIn[d[1]]*u + pr; // rho*u*u + p
    fOut[d[2]] = qIn[d[2]]*u; // rho*v*u
    fOut[d[3]] = qIn[d[3]]*u; // rho*w*u
    fOut[4] = (qIn[4]+pr)*u; // (E+p)*u
  }
}
