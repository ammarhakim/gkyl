// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FiveMomentSrcImpl.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Eigen>

// Makes indexing cleaner
static const unsigned X = 0;
static const unsigned Y = 1;
static const unsigned Z = 2;

static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

#define sq(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))
template <typename T> static T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

static const int COL_PIV_HOUSEHOLDER_QR = 0;
static const int PARTIAL_PIV_LU = 1;

#define fidx(n, c) (3 * (n) + (c))
#define eidx(c) (3 * nFluids + (c))

void
gkylFiveMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcRk3(sd, fd, dt, ff, em);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcTimeCentered(sd, fd, dt, ff, em, staticEm, sigma);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcTimeCenteredDirect2(sd, fd, dt, ff, em, staticEm, sigma);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcTimeCenteredDirect(sd, fd, dt, ff, em, staticEm, sigma);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt,
                       double **ff, double *em, double *staticEm, double *sigma)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcExact(sd, fd, dt, ff, em, staticEm, sigma);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}
