//------------------------------------------------------------------------------
// C++ back-end for inter-species friction for the five-moment (Euler) model.
//------------------------------------------------------------------------------

#include <FiveMomentFrictionSrcImpl.h>
#include <cmath>
#include <Eigen/Eigen>
#include <cassert>
#include <iostream>

static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

#define NN (0)
#define VX (1)
#define VY (2)
#define VZ (3)
#define PP (4)
#define TT (4)

#define sq(x) ((x) * (x))


static inline double
pressure(const double *f, const double gamma)
{
  return (f[ER] - 0.5 * (sq(f[MX])+sq(f[MY])+sq(f[MZ])) / f[RHO]) * (gamma-1);
}


static inline double
energy(const double *f, const double p, const double gamma)
{
  return p / (gamma-1) + 0.5 * (sq(f[MX])+sq(f[MY])+sq(f[MZ])) / f[RHO];
}


static void
calcNu(const FiveMomentFrictionSrcData_t *sd,
       const double * const *fPtrs,
       const double *nuBase,
       double *nu)
{
  const int nFluids = sd->nFluids;

  // Read pre-specified nu_sr for r>s from the flattend 1d array nuBase.
  // FIXME Do this only once.
  int i = 0;
  for (int s=0; s<nFluids; ++s)
  {
    double *nu_s = nu + nFluids * s;
    for (int r=s+1; r<nFluids; ++r)
    {
      nu_s[r] = nuBase[i];
      i += 1;
    }
  }

  // Compute the remaining nu terms that satisfy nu_sr*rho_s == nu_rs*rho_r, and
  // store store the result in the flattened 1d array nu.
  for (int s=0; s<nFluids; ++s)
  {
    const double rho_s = fPtrs[s][RHO];
    double *nu_s = nu + nFluids * s;

    for (int r=0; r<s; ++r)
    {
      const double rho_r = fPtrs[r][RHO];
      const double *nu_r = nu + nFluids * r;
      nu_s[r] = nu_r[s] * rho_r / rho_s;
    }
  }
}


// F = f + dt * L(f), L is an operator
static void
eulerUpdate(const FiveMomentFrictionSrcData_t *sd,
            const FiveMomentFrictionFluidData_t *fd,
            const double dt,
            double *fs,
            double *F,
            const int s,
            const double * const *fPtrsOld,
            const double *nu)
{
  const int nFluids = sd->nFluids;

  const double ms = fd[s].mass;
  const double u_s = fs[MX] / fs[RHO];
  const double v_s = fs[MY] / fs[RHO];
  const double w_s = fs[MZ] / fs[RHO];
  const double T_s = pressure(fs, sd->gasGamma) / fs[RHO] * ms;
  const double *nu_s = nu + nFluids * s;

  double u_s1=u_s, v_s1=v_s, w_s1=w_s, T_s1=T_s;

  for (int r=0; r<nFluids; ++r)
  {
    if (r==s)
      continue;

    const double *fr = fPtrsOld[r];
    const double du = fr[MX] / fr[RHO] - u_s;
    const double dv = fr[MY] / fr[RHO] - v_s;
    const double dw = fr[MZ] / fr[RHO] - w_s;

    u_s1 += dt * nu_s[r] * du;
    v_s1 += dt * nu_s[r] * dv;
    w_s1 += dt * nu_s[r] * dw;

    if (sd->hasPressure)
    {
      const double mr = fd[r].mass;
      const double T_r = pressure(fr, sd->gasGamma) / fr[RHO] * mr;
      const double du2 = sq(du) + sq(dv) + sq(dw);
      const double coeff = dt * nu_s[r] * ms / (ms + mr);
      T_s1 += coeff * ( 2*(T_r-T_s) + (mr/3)*du2 );
    }
  }

  F[RHO] = fs[RHO];
  F[MX] = fs[RHO] * u_s1;
  F[MY] = fs[RHO] * v_s1;
  F[MZ] = fs[RHO] * w_s1;
  F[ER] = energy(F, T_s1 * fs[RHO] / ms, sd->gasGamma);
}


void
gkylFiveMomentFrictionSrcForwardEuler(
  const FiveMomentFrictionSrcData_t *sd,
  const FiveMomentFrictionFluidData_t *fd,
  const double dt,
  double **fPtrs)
{
  const int nFluids = sd->nFluids;

  double old[nFluids*5];
  double *fPtrsOld[nFluids];
  for (int s=0; s<nFluids; ++s)
  {
    double *f = fPtrs[s];
    double *fOld = old + s * 5;
    for (int comp=0; comp<5; ++comp)
      fOld[comp] = f[comp];
    fPtrsOld[s] = fOld;
  }

  double nu[nFluids * nFluids];
  calcNu(sd, fPtrs, sd->nuBase, nu);

  for (int s=0; s<nFluids; ++s)
  {
    double *f = fPtrs[s];
    eulerUpdate(sd, fd, dt, f, f, s, fPtrsOld, nu);
  }
}


void
gkylFiveMomentFrictionSrcTimeCentered(
  const FiveMomentFrictionSrcData_t *sd,
  const FiveMomentFrictionFluidData_t *fd,
  const double dt,
  double **fPtrs)
{
  const int nFluids = sd->nFluids;

  double nu[nFluids * nFluids];
  calcNu(sd, fPtrs, sd->nuBase, nu);

  // Update velocities.
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Constant(nFluids, nFluids, 0.0);
  Eigen::VectorXd rhs_u(nFluids);
  Eigen::VectorXd rhs_v(nFluids);
  Eigen::VectorXd rhs_w(nFluids);

  for (int s=0; s<nFluids; ++s)
  {
    double *fs = fPtrs[s];
    rhs_u(s) = fs[MX] / fs[RHO];
    rhs_v(s) = fs[MY] / fs[RHO];
    rhs_w(s) = fs[MZ] / fs[RHO];

    lhs(s, s) = 1.;
    const double *nu_s = nu + nFluids * s;
    for (int r=0; r<nFluids; ++r)
    {
      if (r==s)
        continue;

      const double half_dt_nu_sr = 0.5 * dt * nu_s[r];
      lhs(s, s) += half_dt_nu_sr;
      lhs(s, r) -= half_dt_nu_sr;
    }
  }

  // Compute velocity at half time-step n+1/2.
#if 0
  Eigen::MatrixXd inv = lhs.inverse();
  Eigen::VectorXd sol_u = inv * rhs_u;
  Eigen::VectorXd sol_v = inv * rhs_v;
  Eigen::VectorXd sol_w = inv * rhs_w;
#else
  Eigen::VectorXd sol_u = lhs.partialPivLu().solve(rhs_u);
  Eigen::VectorXd sol_v = lhs.partialPivLu().solve(rhs_v);
  Eigen::VectorXd sol_w = lhs.partialPivLu().solve(rhs_w);
#endif

  // Update temperature, and immediately, pressure.
  if (sd->hasPressure) {
    lhs.setZero();
    Eigen::VectorXd rhs_T(nFluids);

    double T[nFluids];
    for (int s=0; s<nFluids; ++s)
    {
      const double ms = fd[s].mass;
      const double coeff = 0.5 * dt * ms;

      double *fs = fPtrs[s];
      T[s] = pressure(fs, sd->gasGamma) / fs[RHO] * fd[s].mass;
      rhs_T(s) = T[s];
      lhs(s, s) = 1.;

      const double *nu_s = nu + nFluids * s;
      for (int r=0; r<nFluids; ++r)
      {
        if (r==s)
          continue;

        const double mr = fd[r].mass;
        const double du2 = sq(rhs_u(s)-rhs_u(r)) \
                         + sq(rhs_v(s)-rhs_v(r)) \
                         + sq(rhs_w(s)-rhs_w(r));
        const double coeff_sr = coeff * nu_s[r] / (ms + mr);

        rhs_T(s) += coeff_sr * (mr / 3.) * du2;
        lhs(s, s) += coeff_sr * 2;
        lhs(s, r) -= coeff_sr * 2;
      }
    }

    // Compute pressure at n+1/2 and then at n+1.
    Eigen::VectorXd sol_T = lhs.partialPivLu().solve(rhs_T);
    for (int s=0; s<nFluids; ++s)
    {
      double *f = fPtrs[s];
      f[PP] = (2 * sol_T(s) - T[s]) * f[RHO] / fd[s].mass;
    }
  } else {
    for (int s=0; s<nFluids; ++s)
    {
      double *f = fPtrs[s];
      f[PP] = pressure(f, sd->gasGamma);
    }
  }

  // Compute momentum (from velocity) and total energy at full time-step n+1.
  for (int s=0; s<nFluids; ++s)
  {
    double *f = fPtrs[s];
    f[MX] = 2 * f[RHO] * sol_u(s) - f[MX];
    f[MY] = 2 * f[RHO] * sol_v(s) - f[MY];
    f[MZ] = 2 * f[RHO] * sol_w(s) - f[MZ];
    f[ER] = f[PP]/(sd->gasGamma-1) + 0.5*(sq(f[MX])+sq(f[MY])+sq(f[MZ]))/f[RHO];
  }
}


void
gkylFiveMomentFrictionSrcExact(
  const FiveMomentFrictionSrcData_t *sd,
  const FiveMomentFrictionFluidData_t *fd,
  const double dt,
  double **fPtrs)
{
  assert(sd->nFluids==2);

  const double nu01 = sd->nuBase[0];
  if (nu01==0)
    return;

  double *f0 = fPtrs[0];
  double *f1 = fPtrs[1];

  const double rho0 = f0[RHO];
  const double rho1 = f1[RHO];

  double p0 = pressure(f0, sd->gasGamma);
  double p1 = pressure(f1, sd->gasGamma);

  double du2 = 0;
  {
    const double nu10 = nu01 * rho0 / rho1;
    const double nuSum = nu01 + nu10;
    const double coeff = (std::exp(-nuSum*dt) - 1) / nuSum;
    const double coeff0 = coeff * nu01;
    const double coeff1 = coeff * nu10;

    for(int d=0; d<3; ++d)
    {
      double u0 = f0[MX+d] / rho0;
      double u1 = f1[MX+d] / rho1;
      double du = u1 - u0;
      u0 -= du * coeff0;
      u1 += du * coeff1;
      f0[MX+d] = rho0 * u0;
      f1[MX+d] = rho1 * u1;

      if (sd->hasPressure) {
        du2 += sq(u1-u0);
      }
    }
  }

  if (sd->hasPressure) {
    const double m0 = fd[0].mass, n0 = rho0 / m0;
    const double m1 = fd[1].mass, n1 = rho1 / m1;

    const double a = 2 * nu01 * rho0 / (m0 + m1);
    const double b0 = (m0/6) * du2;
    const double b1 = (m1/6) * du2;

    // Heating in the total pressure.
    double p = p0 + p1;
    p += (b0 + b1) * a * dt;

    const double rate = a * (n0+n1) / (n0*n1);
    const double shift = (n0*b0-n1*b1) / (n0+n1);

    // Temperature difference relaxation.
    double dT = p0/n0 - p1/n1;
    dT = (dT + shift) * std::exp(-rate*dt) - shift;

    p0 = n0 * (p + n1 * dT) / (n0 + n1);
    p1 = p - p0;
  }

  f0[ER] = energy(f0, p0, sd->gasGamma);
  f1[ER] = energy(f1, p1, sd->gasGamma);
}
