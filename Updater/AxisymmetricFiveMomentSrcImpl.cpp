// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for geometric source for axisymmetric five-moment (Euler) model
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <AxisymmetricFiveMomentSrcImpl.h>
#include <cmath>
#include <cassert>
#include <iostream>

static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;


// F = f + dt * L(f), L is an operator
static void
eulerUpdate(const AxisymmetricFiveMomentSrcData_t *sd,
            const double dt,
            const double *xc,
            double *f,
            double *F)
{
  const double radius = xc[0];
  const double gasGamma = sd->gasGamma;

  const double rho = f[RHO];
  const double u = f[MX] / rho;
  const double v = f[MY] / rho;
  const double w = f[MZ] / rho;
  double e_in = sd->hasPressure ? (f[ER] - 0.5 * rho * (u*u + v*v + w*w)) : 0.;

  F[RHO] = f[RHO] - (dt/radius) * (rho*u);
  F[MX]  = f[MX] - (dt/radius) * (rho*u*u - rho*v*v);
  F[MY]  = f[MY] - (dt/radius) * (2*rho*u*v);
  F[MZ]  = f[MZ] - (dt/radius) * (rho*u*w);

  if (sd->hasPressure) {
    e_in -= (dt/radius) * gasGamma * u * e_in;
    F[ER] = e_in + 0.5*(f[MX]*f[MX]+f[MY]*f[MY]+f[MZ]*f[MZ])/f[RHO];
  }

  // if (sd->hasPressure) {
  //   const double E = f[ER];
  //   const double p = (sd->gasGamma - 1) * (E - 0.5 * rho * (u*u + v*v + w*w));
  //   F[ER] = f[ER] - (dt/radius) * (u*(E+p));
  // }
}

inline static void
combine(const int ncomps,
        const double c1,
        double *f1,
        const double c2,
        double *f2,
        double *fOut)
{
  for (int c=0; c<ncomps; ++c)
  {
    fOut[c] = c1 * f1[c] + c2 * f2[c];
  }
}

void
gkylAxisymmetricFiveMomentSrcForwardEuler(
  const AxisymmetricFiveMomentSrcData_t *sd,
  const AxisymmetricFluidData_t *fd,
  const double dt,
  const double *xc,
  double **fPtrs)
{
  for (int s=0; s<sd->nFluids; ++s)
  {
    if (!(fd[s].evolve))
      continue;

    double *f = fPtrs[s];

    eulerUpdate(sd, dt, xc, f, f);
  }
}

void
gkylAxisymmetricFiveMomentSrcRk3(const AxisymmetricFiveMomentSrcData_t *sd,
                                 const AxisymmetricFluidData_t *fd,
                                 const double dt,
                                 const double *xc,
                                 double **fPtrs)
{
  const int ncomps = sd->hasPressure? 5 : 4;
  double f1[ncomps], f2[ncomps];
  for (int s=0; s<sd->nFluids; ++s)
  {
    if (!(fd[s].evolve))
      continue;

    double *f = fPtrs[s];

    // rk stage 1, f1 = F(f) = f + dt*L(f)
    eulerUpdate(sd, dt, xc, f, f1);

    // rk stage 2, f2 = 0.75*f + 0.25 * F(f1)
    eulerUpdate(sd, dt, xc, f1, f1);
    combine(ncomps, 0.75, f, 0.25, f1, f2);

    // rk stage 3, f = (1/3) * f + (2/3) * F(f2)
    eulerUpdate(sd, dt, xc, f2, f2);
    combine(ncomps, 1./3., f, 2./3., f2, f);
  }
}

// F = f + dt * L(f), L is an operator
static void
eulerUpdatePositivity(const AxisymmetricFiveMomentSrcData_t *sd,
            const double dt,
            const double *xc,
            double *f,
            double *F)
{
  const double radius = xc[0];
  const double gasGamma = sd->gasGamma;

  const double rho = f[RHO];
  const double u = f[MX] / rho;
  const double v = f[MY] / rho;
  const double w = f[MZ] / rho;
  double e_in = sd->hasPressure ? (f[ER] - 0.5 * rho * (u*u + v*v + w*w)) : 0.;

  // Explicit update of momentum with rho being a fixed value.
  F[MX]  = f[MX] - (dt/radius) * (rho*u*u - rho*v*v);
  F[MY]  = f[MY] - (dt/radius) * (2*rho*u*v);
  F[MZ]  = f[MZ] - (dt/radius) * (rho*u*w);

  // Exact Update of density and pressure using the new radial velocity as a
  // fixed value.
  const double u1 = f[MX] / rho;
  f[RHO] *= std::exp(-u1*dt/radius);

  if (sd->hasPressure) {
    e_in *= std::exp(-u1*dt/radius * gasGamma);;
    F[ER] = e_in + 0.5*(f[MX]*f[MX]+f[MY]*f[MY]+f[MZ]*f[MZ])/f[RHO];
  }
}

void
gkylAxisymmetricFiveMomentSrcPositivityForwardEuler(
  const AxisymmetricFiveMomentSrcData_t *sd,
  const AxisymmetricFluidData_t *fd,
  const double dt,
  const double *xc,
  double **fPtrs)
{
  for (int s=0; s<sd->nFluids; ++s)
  {
    if (!(fd[s].evolve))
      continue;

    double *f = fPtrs[s];

    eulerUpdatePositivity(sd, dt, xc, f, f);
  }
}
