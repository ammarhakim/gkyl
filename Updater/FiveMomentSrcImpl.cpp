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

static double sq(double x) { return x*x; }

static const int COL_PIV_HOUSEHOLDER_QR = 0;
static const int PARTIAL_PIV_LU = 1;

#define fidx(n, c) (3 * (n) + (c))
#define eidx(c) (3 * nFluids + (c))

void
gkylFiveMomentSrcRk3(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em)
{
  unsigned nFluids = sd->nFluids;
  std::vector<double> f1(5), em1(8), curr(3);
  std::vector<double> f2(5), em2(8);

  // if fluid equations have pressure, then compute initial KE
  double keOld = 0.0;
  if (sd->hasPressure)
  {
    for (unsigned n=0; n<nFluids; ++n)
      keOld = 0.5*(sq(ff[n][MX]) + sq(ff[n][MY]) + sq(ff[n][MZ]))/ff[n][RHO];
  }

  // B fields don't change
  em1[BX] = em[BX]; em1[BY] = em[BY]; em1[BZ] = em[BZ];
  em2[BX] = em[BX]; em2[BY] = em[BY]; em2[BZ] = em[BZ];

  //------------> RK Stage 1
  curr[0] = curr[1] = curr[2] = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  { // update fluids
    double *f = ff[n];
    double qmdt = dt*fd[n].charge/fd[n].mass;
    double qmdte = qmdt/sd->epsilon0;

    f1[MX] = f[MX] + qmdt*(em[EX] + f[MY]*em[BZ] - f[MZ]*em[BY]);
    f1[MY] = f[MY] + qmdt*(em[EY] + f[MZ]*em[BX] - f[MX]*em[BZ]);
    f1[MZ] = f[MZ] + qmdt*(em[EZ] + f[MX]*em[BY] - f[MY]*em[BX]);

    curr[EX] += qmdte*f[MX];
    curr[EY] += qmdte*f[MY];
    curr[EZ] += qmdte*f[MZ];
  }
  // update field
  em1[EX] = em[EX] - curr[EX];
  em1[EY] = em[EY] - curr[EY];
  em1[EZ] = em[EZ] - curr[EZ];

  //------------> RK Stage 2
  curr[0] = curr[1] = curr[2] = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  { // update fluids
    double *f = ff[n];
    double qmdt = dt*fd[n].charge/fd[n].mass;
    double qmdte = qmdt/sd->epsilon0;

    f2[MX] = 0.75*f[MX] + 0.25*(f1[MX] + qmdt*(em1[EX] + f1[MY]*em1[BZ] - f1[MZ]*em1[BY]));
    f2[MY] = 0.75*f[MY] + 0.25*(f1[MY] + qmdt*(em1[EY] + f1[MZ]*em1[BX] - f1[MX]*em1[BZ]));
    f2[MZ] = 0.75*f[MY] + 0.25*(f1[MZ] + qmdt*(em1[EZ] + f1[MX]*em1[BY] - f1[MY]*em1[BX]));

    curr[EX] += qmdte*f1[MX];
    curr[EY] += qmdte*f1[MY];
    curr[EZ] += qmdte*f1[MZ];
  }
  // update field
  em2[EX] = 0.75*em[EX] + 0.25*(em1[EX] - curr[EX]);
  em2[EY] = 0.75*em[EY] + 0.25*(em1[EY] - curr[EY]);
  em2[EZ] = 0.75*em[EZ] + 0.25*(em1[EZ] - curr[EZ]);

  //------------> RK Stage 3
  double one3 = 1.0/3.0, two3 = 2.0/3.0;
  curr[0] = curr[1] = curr[2] = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  { // update fluids
    double *f = ff[n];
    double qmdt = dt*fd[n].charge/fd[n].mass;
    double qmdte = qmdt/sd->epsilon0;

    f[MX] = one3*f[MX] + two3*(f2[MX] + qmdt*(em2[EX] + f2[MY]*em2[BZ] - f2[MZ]*em2[BY]));
    f[MY] = one3*f[MY] + two3*(f2[MY] + qmdt*(em2[EY] + f2[MZ]*em2[BX] - f2[MX]*em2[BZ]));
    f[MZ] = one3*f[MY] + two3*(f2[MZ] + qmdt*(em2[EZ] + f2[MX]*em2[BY] - f2[MY]*em2[BX]));

    curr[EX] += qmdte*f2[MX];
    curr[EY] += qmdte*f2[MY];
    curr[EZ] += qmdte*f2[MZ];
  }
  // update field
  em[EX] = one3*em[EX] + two3*(em1[EX] - curr[EX]);
  em[EY] = one3*em[EY] + two3*(em1[EY] - curr[EY]);
  em[EZ] = one3*em[EZ] + two3*(em1[EZ] - curr[EZ]);

  // if fluid equations have pressure, update energy
  if (sd->hasPressure)
  {
    for (unsigned n=0; n<nFluids; ++n)
    {
      double keNew = 0.5*(sq(ff[n][MX]) + sq(ff[n][MY]) + sq(ff[n][MZ]))/ff[n][RHO];
      ff[n][ER] += keNew-keOld;
    }
  }  
}

void
gkylFiveMomentSrcTimeCentered(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  unsigned nFluids = sd->nFluids;
  double dt1 = 0.5 * dt;
  double dt2 = 0.5 * dt / sd->epsilon0;
  
  double zeros[6] = {0.};
  if (!(sd->hasStatic))
  {
    staticEm = zeros;
  }

  Eigen::MatrixXd lhs = Eigen::MatrixXd::Constant(3*nFluids+3, 3*nFluids+3, 0.0);
  Eigen::VectorXd rhs(3*nFluids+3);

  //------------> fill elements corresponding to each fluid
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;
    double qbym2 = sq(qbym);

    // eqn. for X-component of current
    lhs(fidx(n,X), fidx(n,X)) = 1.0;
    lhs(fidx(n,X), fidx(n,Y)) = -dt1*qbym*(em[BZ]+staticEm[BZ]);
    lhs(fidx(n,X), fidx(n,Z)) = dt1*qbym*(em[BY]+staticEm[BY]);
    lhs(fidx(n,X), eidx(X)) = -dt1*qbym2*f[RHO];

    // eqn. for Y-component of current
    lhs(fidx(n,Y), fidx(n,X)) = dt1*qbym*(em[BZ]+staticEm[BZ]);
    lhs(fidx(n,Y), fidx(n,Y)) = 1.0;
    lhs(fidx(n,Y), fidx(n,Z)) = -dt1*qbym*(em[BX]+staticEm[BX]);
    lhs(fidx(n,Y), eidx(Y)) = -dt1*qbym2*f[RHO];

    // eqn. for Z-component of current
    lhs(fidx(n,Z), fidx(n,X)) = -dt1*qbym*(em[BY]+staticEm[BY]);
    lhs(fidx(n,Z), fidx(n,Y)) = dt1*qbym*(em[BX]+staticEm[BX]);
    lhs(fidx(n,Z), fidx(n,Z)) = 1.0;
    lhs(fidx(n,Z), eidx(Z)) = -dt1*qbym2*f[RHO];

    // fill corresponding RHS elements
    rhs(fidx(n,X)) = qbym*f[MX];
    rhs(fidx(n,Y)) = qbym*f[MY];
    rhs(fidx(n,Z)) = qbym*f[MZ];

    // add gravity source term for current
    rhs(fidx(n, sd->gravityDir)) += qbym*f[RHO]*sd->gravity*dt1;

    // set current contribution to electric field equation
    lhs(eidx(X), fidx(n,X)) = dt2;
    lhs(eidx(Y), fidx(n,Y)) = dt2;
    lhs(eidx(Z), fidx(n,Z)) = dt2;
  }

  // fill in elements for electric field equations
  Eigen::VectorXd sol;
  lhs(eidx(EX), eidx(EX)) = 1.0;
  lhs(eidx(EY), eidx(EY)) = 1.0;
  lhs(eidx(EZ), eidx(EZ)) = 1.0;

  rhs(eidx(EX)) = em[EX];
  rhs(eidx(EY)) = em[EY];
  rhs(eidx(EZ)) = em[EZ];

  // invert to find solution
  if (sd->linSolType == COL_PIV_HOUSEHOLDER_QR)
    sol = lhs.colPivHouseholderQr().solve(rhs);
  else if (sd->linSolType == PARTIAL_PIV_LU)
    sol = lhs.partialPivLu().solve(rhs);
  else
    { /* can not happen */ }

  //------------> update solution for fluids
  double chargeDens = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;

    chargeDens += qbym*f[RHO];

    double keOld = 0.5*(sq(f[MX]) + sq(f[MY]) + sq(f[MZ]))/f[RHO];
    f[MX] = 2*sol(fidx(n,X))/qbym - f[MX];
    f[MY] = 2*sol(fidx(n,Y))/qbym - f[MY];
    f[MZ] = 2*sol(fidx(n,Z))/qbym - f[MZ];
    if (sd->hasPressure)
    {
      double keNew = 0.5*(sq(f[MX]) + sq(f[MY]) + sq(f[MZ]))/f[RHO];
      f[ER] += keNew-keOld;
    }
  }

  //------------> update electric field
  em[EX] = 2*sol(eidx(X)) - em[EX];
  em[EY] = 2*sol(eidx(Y)) - em[EY];
  em[EZ] = 2*sol(eidx(Z)) - em[EZ];

  //------------> update correction potential
  double crhoc = sd->chi_e*chargeDens/sd->epsilon0;
  em[PHIE] += dt*crhoc;
}
