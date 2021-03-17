// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for geometric source for axisymmetric perfectly-hyperbolic
// Maxwell's equations
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <AxisymmetricPhMaxwellSrcImpl.h>

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

// F = f + dt * L(f), L is an operator
static void
eulerUpdate(const AxisymmetricPhMaxwellSrcData_t *sd,
            const double dt,
            const double *xc,
            double *f,
            double *F)
{
  const double radius = xc[0];
  const double epsilon0 = sd->epsilon0;
  const double mu0 = sd->mu0;
  const double c2 = 1. / (epsilon0 * mu0);
  const double chi_e = sd->chi_e;
  const double chi_m = sd->chi_m;

  F[EZ] = f[EZ] + (dt/radius) * c2 * f[BY];
  F[BZ] = f[BZ] - (dt/radius) * f[EY];
  F[PHIE] = f[PHIE] - (dt/radius) * chi_e * f[EX];
  F[PHIM] = f[PHIM] - (dt/radius) * chi_m * c2 * f[BX];

  if (!(f==F))
  {
    F[EX] = f[EX];
    F[EY] = f[EY];
    F[BX] = f[BX];
    F[BY] = f[BY];
  }
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
gkylAxisymmetricPhMaxwellSrcForwardEuler(
  const AxisymmetricPhMaxwellSrcData_t *sd,
  const double dt,
  const double *xc,
  double *f)
{
  eulerUpdate(sd, dt, xc, f, f);
}

void
gkylAxisymmetricPhMaxwellSrcRk3(const AxisymmetricPhMaxwellSrcData_t *sd,
                                const double dt,
                                const double *xc,
                                double *f)
{
  const int ncomps = 8;
  double f1[ncomps], f2[ncomps];

  // rk stage 1, f1 = F(f) = f + dt*L(f)
  eulerUpdate(sd, dt, xc, f, f1);

  // rk stage 2, f2 = 0.75*f + 0.25 * F(f1)
  eulerUpdate(sd, dt, xc, f1, f1);
  combine(ncomps, 0.75, f, 0.25, f1, f2);

  // rk stage 3, f = (1/3) * f + (2/3) * F(f2)
  eulerUpdate(sd, dt, xc, f2, f2);
  combine(ncomps, 1./3., f, 2./3., f2, f);
}

