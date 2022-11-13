//------------------------------------------------------------------------------
// C++ back-end for user-supplied source terms for the five-moment (Euler) model.
//------------------------------------------------------------------------------

#include <FiveMomentUserSrcImpl.h>

#define RHO (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)

#define NN (0)
#define TT (1)

#define sq(x) ((x) * (x))

// q1 = q + dt * L(q), L is an operator
// sd: solver parameters
// fd: fluid species parameters
// dt: time step size
// q: input fluid species data
// q1: output fluid species data after update
// S: source data
static void
eulerUpdate(const FiveMomentUserSrcData_t *sd,
            const FiveMomentUserSrcFluidData_t *fd,
            const double dt,
            double *q,
            double *q1,
            const double *S)
{
  double m = fd->mass, gamma = sd->gasGamma, kB = sd->kBoltzmann;

  double rho0 = q[RHO];
  double n0 = rho0 / m;
  double u = q[MX] / rho0;
  double v = q[MY] / rho0;
  double w = q[MZ] / rho0;
  double v2 = u*u + v*v + w*w;
  double T0 = (q[ER] - 0.5 * rho0 * v2) * (gamma-1.0) / n0 / kB;

  double n1 = n0 + dt * S[NN];
  double T1 = T0 + dt * S[TT];

  double rho1 = n0 * m;
  q[RHO] = rho1;
  q1[MX] = rho1 * u;
  q1[MY] = rho1 * v;
  q1[MZ] = rho1 * w;
  q1[ER] = n1 * kB * T1 / (gamma-1.0) + 0.5 * rho1 * v2;
}

void
gkylFiveMomentUserSrcForwardEuler(
  const FiveMomentUserSrcData_t *sd,
  const FiveMomentUserSrcFluidData_t *fd,
  const double dt,
  double *qPtrs[],
  const double *sourcePtrs[])
{
  for (int s=0; s<sd->nFluids; ++s)
  {
    double *q = qPtrs[s];
    const double *S = sourcePtrs[s];
    eulerUpdate(sd, fd+s, dt, q, q, S);
  }
}
