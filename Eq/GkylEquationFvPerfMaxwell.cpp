// Gkyl ------------------------------------------------------------------------
//
// Perfectly-Hyperbolic Maxwell's equations object for wave-propagation
// finite-volume method
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylEquationFvPerfMaxwell.h>

__constant__ const int dirShufflePerfMaxwell[3][7] = {
  {0,1,2,3,4,5,6},
  {0,2,3,1,5,6,4},
  {0,3,1,2,6,4,5}
};

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

  double delta[8];
#pragma unroll
  for (int c=0; c<8; c++) {
    delta[c] = qr[c]-ql[c];
  }

  // compute projections of jump (generated from Maxima)
  const double a1 = 0.5*(delta[d[3]]-delta[7]*c1);
  const double a2 = 0.5*(delta[7]*c1+delta[d[3]]);
  const double a3 = 0.5*(delta[d[0]]-delta[6]*c);
  const double a4 = 0.5*(delta[6]*c+delta[d[0]]);
  const double a5 = 0.5*(delta[d[1]]-delta[d[5]]*c);
  const double a6 = 0.5*(delta[d[4]]*c+delta[d[2]]);
  const double a7 = 0.5*(delta[d[5]]*c+delta[d[1]]);
  const double a8 = 0.5*(delta[d[2]]-delta[d[4]]*c);

  // set waves to 0.0 as most entries vanish
#pragma unroll
  for (int i=0; i<8*6; i++) {
    waves[i] = 0;
  }
  // memset(waves, 0, sizeof(double)*48);

  // wave 1:
  double *w = waves;
  w[d[3]] = a1;
  w[7] = -a1*c;
  speeds[0] = -c*cb;

  // wave 2:
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
  w[d[4]] = a6*c1;
  w[d[5]] = -a5*c1;
  speeds[4] = -c;

  // wave 6: (two waves with EV c, c lumped into one)
  w += 8;
  w[d[1]] = a7;
  w[d[2]] = a8;
  w[d[4]] = -a8*c1;
  w[d[5]] = a7*c1;
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
  {
    const double *w1p = waves+1*8;
    const double *w2p = waves+3*8;
    const double *w3p = waves+5*8;
    const double s1p = speeds[1];
    const double s2p = speeds[3];
    const double s3p = speeds[5];

    for (int i=0; i<8; i++)
    {
      apdq[i] = s1p*w1p[i] + s2p*w2p[i] + s3p*w3p[i];
    }
  }

  {
    const double *w1m = waves+0*8;
    const double *w2m = waves+2*8;
    const double *w3m = waves+4*8;
    const double s1m = speeds[0];
    const double s2m = speeds[2];
    const double s3m = speeds[4];

    for (int i=0; i<8; i++)
    {
      amdq[i] = s1m*w1m[i] + s2m*w2m[i] + s3m*w3m[i];
    }
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
