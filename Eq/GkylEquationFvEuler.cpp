// Gkyl ------------------------------------------------------------------------
//
// Euler equation object for wave-propagation finite-volume method
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylEquationFvEuler.h>

__host__ __device__ inline double Euler_pressure(
    const double gasGamma,
    const double * __restrict__ q)
{
  return (gasGamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
}

__host__ __device__ void Euler_rp(
    const void * __restrict__ eqn,
    const int dir,
    const double * __restrict__ ql,
    const double * __restrict__ qr,
    double * __restrict__ waves,
    double * __restrict__ speeds)
{
  const GkylEquationFvEuler_t *eulerEqn = (GkylEquationFvEuler_t *)eqn;
  const int *d = dirShuffle[dir];
  const double gasGamma = eulerEqn->gasGamma;

  double *wv = &waves[0];
#pragma unroll
  for (int c=0; c<5; c++) {
    wv[c] = qr[c]-ql[c];
  }

  double pr = Euler_pressure(gasGamma, ql);
  double u = ql[d[1]]/ql[0];
  double cs = sqrt(gasGamma*pr/ql[0]);
  double sl = u+cs;

  pr = Euler_pressure(gasGamma, qr);
  u = qr[d[1]]/qr[0];
  cs = sqrt(gasGamma*pr/qr[0]);
  double sr = u+cs   ;

  speeds[0] = 0.5*(sl+sr);
}
__device__ static const rp_t p_Euler_rp = &Euler_rp;

__host__ __device__ void Euler_qFluctuations(
    const void * __restrict__ eqn,
    const int dir,
    const double * __restrict__ ql,
    const double * __restrict__ qr,
    const double * __restrict__ waves,
    const double * __restrict__ speeds,
    double * __restrict__ amdq,
    double * __restrict__ apdq) {
  const GkylEquationFvEuler_t *eulerEqn = (GkylEquationFvEuler_t *)eqn;
  const int *d = dirShuffle[dir];
  const double gasGamma = eulerEqn->gasGamma;

  // XXX TODO set size to numEquations
  double fl[10];
  double fr[10];

  double pr = Euler_pressure(gasGamma, ql);
  double u = ql[d[1]]/ql[0];
  fl[0] = ql[d[1]]; // rho*u
  fl[d[1]] = ql[d[1]]*u + pr; // rho*u*u + p
  fl[d[2]] = ql[d[2]]*u; // rho*v*u
  fl[d[3]] = ql[d[3]]*u; // rho*w*u
  fl[4] = (ql[4]+pr)*u; // (E+p)*u
  double absMaxsl = abs(u)+sqrt(gasGamma*pr/ql[0]);

  pr = Euler_pressure(gasGamma, qr);
  u = qr[d[1]]/qr[0];
  fr[0] = qr[d[1]]; // rho*u
  fr[d[1]] = qr[d[1]]*u + pr; // rho*u*u + p
  fr[d[2]] = qr[d[2]]*u; // rho*v*u
  fr[d[3]] = qr[d[3]]*u; // rho*w*u
  fr[4] = (qr[4]+pr)*u; // (E+p)*u
  double absMaxsr = abs(u)+sqrt(gasGamma*pr/qr[0]);

  double absMaxs = max(absMaxsl, absMaxsr);

  amdq[0] = 0.5*(fr[0]-fl[0] - absMaxs*(qr[0]-ql[0]));
  amdq[d[1]] = 0.5*(fr[d[1]]-fl[d[1]] - absMaxs*(qr[d[1]]-ql[d[1]]));
  amdq[d[2]] = 0.5*(fr[d[2]]-fl[d[2]] - absMaxs*(qr[d[2]]-ql[d[2]]));
  amdq[d[3]] = 0.5*(fr[d[3]]-fl[d[3]] - absMaxs*(qr[d[3]]-ql[d[3]]));
  amdq[4] = 0.5*(fr[4]-fl[4] - absMaxs*(qr[4]-ql[4]));

  apdq[0]  = 0.5*(fr[0]-fl[0] + absMaxs*(qr[0]-ql[0]));
  apdq[d[1]] = 0.5*(fr[d[1]]-fl[d[1]] + absMaxs*(qr[d[1]]-ql[d[1]]));
  apdq[d[2]] = 0.5*(fr[d[2]]-fl[d[2]] + absMaxs*(qr[d[2]]-ql[d[2]]));
  apdq[d[3]] = 0.5*(fr[d[3]]-fl[d[3]] + absMaxs*(qr[d[3]]-ql[d[3]]));
  apdq[4] = 0.5*(fr[4]-fl[4] + absMaxs*(qr[4]-ql[4]));
}
__device__ static const qFluctuations_t p_Euler_qFluctuations = 
  &Euler_qFluctuations;

__host__ __device__ void Euler_flux(
    const void * __restrict__ eqn,
    const int dir,
    const double * __restrict__ qIn,
    double * __restrict__ fOut)
{
  const GkylEquationFvEuler_t *eulerEqn = (GkylEquationFvEuler_t *)eqn;
  const int *d = dirShuffle[dir];
  double pr = Euler_pressure(eulerEqn->gasGamma, qIn);
  double u = qIn[d[1]]/qIn[0];
  fOut[0] = qIn[d[1]]; // rho*u.
  fOut[d[1]] = qIn[d[1]]*u + pr; // rho*u*u + p
  fOut[d[2]] = qIn[d[2]]*u; // rho*v*u
  fOut[d[3]] = qIn[d[3]]*u; // rho*w*u
  fOut[4] = (qIn[4]+pr)*u; // (E+p)*u
}
__device__ static const flux_t p_Euler_flux = &Euler_flux;

GkylEquationFvEuler_t *new_EquationFvEulerOnHost(const double gasGamma)
{
  const int numEquations = 5; 
  const int numWaves = 1;

  GkylEquationFvEuler_t *eulerEqn = new GkylEquationFvEuler_t;

  eulerEqn->numWaves = numWaves;
  eulerEqn->numEquations = numEquations;
  eulerEqn->gasGamma = gasGamma;

  return eulerEqn;
}

GkylEquationFv_t *new_EquationFvEulerOnDevice(const double gasGamma)
{
  const int numEquations = 5; 
  const int numWaves = 1;

  GkylEquationFv_t *eqn = new GkylEquationFv_t;
  GkylEquationFvEuler_t *eulerEqn = new GkylEquationFvEuler_t;

  eulerEqn->numWaves = numWaves;
  eulerEqn->numEquations = numEquations;
  eulerEqn->gasGamma = gasGamma;

  eqn->numWaves = numWaves;
  eqn->numEquations = numEquations;
  eqn->equation = Gkyl::CudaUtils<GkylEquationFvEuler_t>::allocAndCopyToDevice(
      eulerEqn);
  
  cudacall(cudaMemcpyFromSymbol(&eqn->equationRp, p_Euler_rp, sizeof(rp_t)));
  cudacall(cudaMemcpyFromSymbol(
        &eqn->equationQFluctuations, p_Euler_qFluctuations,
        sizeof(qFluctuations_t)));
  cudacall(cudaMemcpyFromSymbol(
        &eqn->equationFlux, p_Euler_flux, sizeof(flux_t)));

  return Gkyl::CudaUtils<GkylEquationFv_t>::allocAndCopyToDevice(eqn);
}
