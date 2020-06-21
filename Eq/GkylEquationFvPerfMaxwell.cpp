// Gkyl ------------------------------------------------------------------------
//
// Perfectly-Hyperbolic Maxwell's equations object for wave-propagation
// finite-volume method
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylEquationFvPerfMaxwell.h>

__host__ __device__ void PerfMaxwell_rp(
    const void * __restrict__ eqn,
    const int dir,
    const double * __restrict__ ql,
    const double * __restrict__ qr,
    double * __restrict__ waves,
    double * __restrict__ speeds)
{
  const GkylEquationFvPerfMaxwell_t *perfMaxwellEqn = (GkylEquationFvPerfMaxwell_t *)eqn;
  const int *d = dirShufflePerfMaxwell[dir];
  const double c = perfMaxwellEqn->lightSpeed;
  const double ce = perfMaxwellEqn->elcErrorSpeedFactor;
  const double cb = perfMaxwellEqn->mgnErrorSpeedFactor;

  const double c1 = 1/c;
  double v = c;
  double a = 0.5 * (v/c-1);

  double delta[8];
#pragma unroll
  for (int c=0; c<8; c++) {
    delta[c] = qr[c]-ql[c];
  }

  // compute projections of jump (generated from Maxima)
  double a1 = 0.5*(delta[d[3]]-delta[7]*c1);
  double a2 = 0.5*(delta[7]*c1+delta[d[3]]);
  double a3 = 0.5*(delta[d[0]]-delta[6]*c);
  double a4 = 0.5*(delta[6]*c+delta[d[0]]);
  double a5 = 0.5*(delta[d[1]]-delta[d[5]]*c);
  double a6 = 0.5*(delta[d[4]]*c+delta[d[2]]);
  double a7 = 0.5*(delta[d[5]]*c+delta[d[1]]);
  double a8 = 0.5*(delta[d[2]]-delta[d[4]]*c);

  // set waves to 0.0 as most entries vanish
#pragma unroll
  for (int i=0; i<8*6; i++) {
    waves[i] = 0;
  }

  // wave 1:
  double *w = waves;
  w[d[3]] = a1;
  w[7] = -a1*c;
  speeds[0] = -c*cb;

  // wave 2:
  // XXX jump
  w += 8;
  w[d[3]] = a2;
  w[7] = a2*c;
  speeds[1] = c*cb;

  // wave 3:
  w += 8;
  w[d[0]] = a3;
  w[6] = -a3*c1;
  speeds[2] = -c*ce;

  // wave 4:;
  w += 8;
  w[d[0]] = a4;
  w[6] = a4*c1;
  speeds[3] = c*ce;

  // wave 5: (two waves with EV -c, -c lumped into one)
  w += 8;
  w[d[1]] = a5;
  w[d[2]] = a6;
  w[d[4]] = a6*c1 + a * delta[d[4]];
  w[d[5]] = -a5*c1 + a * delta[d[5]];
  speeds[4] = -c;

  // wave 6: (two waves with EV c, c lumped into one)
  w += 8;
  w[d[1]] = a7;
  w[d[2]] = a8;
  w[d[4]] = -a8*c1 + a * delta[d[4]];
  w[d[5]] = a7*c1 + a * delta[d[5]];
  speeds[5] = c;
}
__device__ static const rp_t p_PerfMaxwell_rp = &PerfMaxwell_rp;

__host__ __device__ void PerfMaxwell_qFluctuations(
    const void * __restrict__ eqn,
    const int dir,
    const double * __restrict__ ql,
    const double * __restrict__ qr,
    const double * __restrict__ waves,
    const double * __restrict__ speeds,
    double * __restrict__ amdq,
    double * __restrict__ apdq)
{
  // FIXME check codes
  double c = speeds[6];
  double v = c;
  const int *d = dirShufflePerfMaxwell[dir]; // shuffle indices for `dir`
  // for _,i in ipairs({2, 3, 5, 6}) do
  // for _,i in ipairs({5, 6}) do
  for (int i=4; i<6; i++) {
     apdq[d[i]] = apdq[d[i]] + 0.5 * (v-c) * (qr[d[i]] - ql[d[i]]);
     amdq[d[i]] = amdq[d[i]] - 0.5 * (v-c) * (qr[d[i]] - ql[d[i]]);
  }

  v = c;
  d = dirShufflePerfMaxwell[dir]; // shuffle indices for `dir`
  // for _,i in ipairs({2, 3, 5, 6}) do
  // for _,i in ipairs({5, 6}) do
  // for _,i in ipairs({2, 3}) do
  for (int i=1; i<3; i++) {
     apdq[d[i]] = apdq[d[i]] + 0.5 * (v-c) * (qr[d[i]] - ql[d[i]]);
     amdq[d[i]] = amdq[d[i]] - 0.5 * (v-c) * (qr[d[i]] - ql[d[i]]);
  }
}
__device__ static const qFluctuations_t p_PerfMaxwell_qFluctuations = 
  &PerfMaxwell_qFluctuations;

__host__ __device__ void PerfMaxwell_flux(
    const void * __restrict__ eqn,
    const int dir,
    const double * __restrict__ qIn,
    double * __restrict__ fOut)
{
  const GkylEquationFvPerfMaxwell_t *perfMaxwellEqn = (GkylEquationFvPerfMaxwell_t *)eqn;
  const int *d = dirShufflePerfMaxwell[dir];
  const double c = perfMaxwellEqn->lightSpeed;
  const double ce = perfMaxwellEqn->elcErrorSpeedFactor;
  const double cb = perfMaxwellEqn->mgnErrorSpeedFactor;
  double c2 = c*c;

  fOut[d[0]] = ce*c2*qIn[6];
  fOut[d[1]] = c2*qIn[d[5]];
  fOut[d[2]] = -c2*qIn[d[4]];
  fOut[d[3]] = cb*qIn[7];
  fOut[d[4]] = -qIn[d[2]];
  fOut[d[5]] = qIn[d[1]];
  fOut[6] = ce*qIn[d[0]];
  fOut[7] = cb*c2*qIn[d[3]];
}
__device__ static const flux_t p_PerfMaxwell_flux = &PerfMaxwell_flux;

GkylEquationFvPerfMaxwell_t *new_EquationFvPerfMaxwellOnHost(
    const double lightSpeed, const double elcErrorSpeedFactor,
    const double mgnErrorSpeedFactor)
{
  const int numEquations = 8; 
  const int numWaves = 6;

  GkylEquationFvPerfMaxwell_t *perfMaxwellEqn = new GkylEquationFvPerfMaxwell_t;

  perfMaxwellEqn->numWaves = numWaves;
  perfMaxwellEqn->numEquations = numEquations;
  perfMaxwellEqn->lightSpeed = lightSpeed;
  perfMaxwellEqn->elcErrorSpeedFactor = elcErrorSpeedFactor;
  perfMaxwellEqn->mgnErrorSpeedFactor = mgnErrorSpeedFactor;

  return perfMaxwellEqn;
}

GkylEquationFv_t *new_EquationFvPerfMaxwellOnDevice(
    const double lightSpeed, const double elcErrorSpeedFactor,
    const double mgnErrorSpeedFactor)
{
  const int numEquations = 8; 
  const int numWaves = 6;

  GkylEquationFv_t *eqn = new GkylEquationFv_t;
  GkylEquationFvPerfMaxwell_t *perfMaxwellEqn = new GkylEquationFvPerfMaxwell_t;

  perfMaxwellEqn->numWaves = numWaves;
  perfMaxwellEqn->numEquations = numEquations;
  perfMaxwellEqn->lightSpeed = lightSpeed;
  perfMaxwellEqn->elcErrorSpeedFactor = elcErrorSpeedFactor;
  perfMaxwellEqn->mgnErrorSpeedFactor = mgnErrorSpeedFactor;

  eqn->numWaves = numWaves;
  eqn->numEquations = numEquations;
  eqn->equation = Gkyl::CudaUtils<GkylEquationFvPerfMaxwell_t>::allocAndCopyToDevice(
      perfMaxwellEqn);
  
  cudacall(cudaMemcpyFromSymbol(&eqn->equationRp, p_PerfMaxwell_rp, sizeof(rp_t)));
  cudacall(cudaMemcpyFromSymbol(
        &eqn->equationQFluctuations, p_PerfMaxwell_qFluctuations,
        sizeof(qFluctuations_t)));
  cudacall(cudaMemcpyFromSymbol(
        &eqn->equationFlux, p_PerfMaxwell_flux, sizeof(flux_t)));

  return Gkyl::CudaUtils<GkylEquationFv_t>::allocAndCopyToDevice(eqn);
}
