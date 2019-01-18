// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <TenMomentSrcImpl.h>
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

static const unsigned P11 = 4;
static const unsigned P12 = 5;
static const unsigned P13 = 6;
static const unsigned P22 = 7;
static const unsigned P23 = 8;
static const unsigned P33 = 9;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

static const int COL_PIV_HOUSEHOLDER_QR = 0;
static const int PARTIAL_PIV_LU = 1;

void
gkylTenMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em)
{
}

static std::vector<double>
gkylPressureTensorSrcTimeCenteredPre(
    MomentSrcData_t *sd, FluidData_t *fd, const double dt, double **ff,
    double *em, double *staticEm)
{
  double Bx = em[BX];
  double By = em[BY];
  double Bz = em[BZ];
  if (sd->hasStatic)
  {
    Bx += staticEm[BX];
    By += staticEm[BY];
    Bz += staticEm[BZ];
  }

  unsigned nFluids = sd->nFluids;
  Eigen::MatrixXd prLhs = Eigen::MatrixXd::Constant(6, 6, 0.0);
  Eigen::VectorXd prRhs(6);
  std::vector<double> prTen(6 * nFluids);
  double dt1 = 0.5 * dt;
  // compute increments from pressure tensor terms
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];

    // assemble LHS terms
    prLhs.noalias() = Eigen::MatrixXd::Identity(6, 6);
    if (fd[n].evolve)
    {
      double qbym = fd[n].charge / fd[n].mass;

      prLhs(0,1) = -2*dt1*qbym*Bz;
      prLhs(0,2) = 2*dt1*qbym*By;
      prLhs(1,0) = dt1*qbym*Bz;
      prLhs(1,2) = -dt1*qbym*Bx;
      prLhs(1,3) = -dt1*qbym*Bz;
      prLhs(1,4) = dt1*qbym*By;
      prLhs(2,0) = -dt1*qbym*By;
      prLhs(2,1) = dt1*qbym*Bx;
      prLhs(2,4) = -dt1*qbym*Bz;
      prLhs(2,5) = dt1*qbym*By;
      prLhs(3,1) = 2*dt1*qbym*Bz;
      prLhs(3,4) = -2*dt1*qbym*Bx;
      prLhs(4,1) = -dt1*qbym*By;
      prLhs(4,2) = dt1*qbym*Bz;
      prLhs(4,3) = dt1*qbym*Bx;
      prLhs(4,5) = -dt1*qbym*Bx;
      prLhs(5,2) = -2*dt1*qbym*By;
      prLhs(5,4) = 2*dt1*qbym*Bx;
    }

    // RHS matrix
    prRhs[0] = f[P11] - f[MX]*f[MX]/f[RHO];
    prRhs[1] = f[P12] - f[MX]*f[MY]/f[RHO];
    prRhs[2] = f[P13] - f[MX]*f[MZ]/f[RHO];
    prRhs[3] = f[P22] - f[MY]*f[MY]/f[RHO];
    prRhs[4] = f[P23] - f[MY]*f[MZ]/f[RHO];
    prRhs[5] = f[P33] - f[MZ]*f[MZ]/f[RHO];

    // solve to compute increment in pressure tensor
    Eigen::VectorXd prSol;
    if (sd->linSolType == COL_PIV_HOUSEHOLDER_QR)
      prSol = prLhs.colPivHouseholderQr().solve(prRhs);
    else if (sd->linSolType == PARTIAL_PIV_LU)
      prSol = prLhs.partialPivLu().solve(prRhs);
    else
    { /* can not happen */ }

    // update solution
    for (unsigned i=0; i<6; ++i)
    {
      prTen[6*n+i] = 2*prSol[i]-prRhs[i];
    }
  }

  return prTen;
}

static std::vector<double>
gkylPressureTensorSrcTimeCenteredDirectPre(
    MomentSrcData_t *sd, FluidData_t *fd, const double dt, double **ff,
    double *em, double *staticEm)
{
  double Bx = em[BX];
  double By = em[BY];
  double Bz = em[BZ];
  if (sd->hasStatic)
  {
    Bx += staticEm[BX];
    By += staticEm[BY];
    Bz += staticEm[BZ];
  }

  unsigned nFluids = sd->nFluids;
  std::vector<double> prTen(6 * nFluids);
  std::vector<double> prRhs(6);
  std::vector<double> prSol(6);
  double dt1 = 0.5 * dt;
  double dtsq = dt1 * dt1;
  double dt3 = dtsq * dt1;
  double dt4 = dt3 * dt1;
  double Bx2 = Bx * Bx;
  double Bx3 = Bx * Bx2;
  double Bx4 = Bx * Bx3;
  double By2 = By * By;
  double By3 = By * By2;
  double By4 = By * By3;
  double Bz2 = Bz * Bz;
  double Bz3 = Bz * Bz2;
  double Bz4 = Bz * Bz3;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    if (!fd[n].evolve)
      continue;
    double *f = ff[n];
    double qbym = fd[n].charge / fd[n].mass;

    prRhs[0] = f[P11] - f[MX] * f[MX] / f[RHO];
    prRhs[1] = f[P12] - f[MX] * f[MY] / f[RHO];
    prRhs[2] = f[P13] - f[MX] * f[MZ] / f[RHO];
    prRhs[3] = f[P22] - f[MY] * f[MY] / f[RHO];
    prRhs[4] = f[P23] - f[MY] * f[MZ] / f[RHO];
    prRhs[5] = f[P33] - f[MZ] * f[MZ] / f[RHO];

    double qb2 = qbym * qbym;
    double qb3 = qbym * qb2;
    double qb4 = qbym * qb3;
    double d = 1 + 5*(Bx2 + By2 + Bz2)*dtsq*qb2 + 4*(Bx2 + By2 + Bz2)*(Bx2 + By2 + Bz2)*dt4*qb4;

    prSol[0] = (prRhs[0] + 2*dt1*(Bz*prRhs[1] - By*prRhs[2])*qbym + dtsq*(5*Bx2*prRhs[0] + 2*Bx*(By*prRhs[1] + Bz*prRhs[2]) + Bz2*(3*prRhs[0] + 2*prRhs[3]) - 4*By*Bz*prRhs[4] + By2*(3*prRhs[0] + 2*prRhs[5]))*qb2 + 2*dt3*(4*Bx2*(Bz*prRhs[1] - By*prRhs[2]) - (By2 + Bz2)*(-(Bz*prRhs[1]) + By*prRhs[2]) - 3*Bx*(By2*prRhs[4] - Bz2*prRhs[4] + By*Bz*(-prRhs[3] + prRhs[5])))*qb3 + 2*dt4*(2*Bx4*prRhs[0] + 4*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 + Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + (By2 + Bz2)*(Bz2*(prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[5])) + Bx2*(4*By*Bz*prRhs[4] + By2*(3*prRhs[3] + prRhs[5]) + Bz2*(prRhs[3] + 3*prRhs[5])))*qb4)/d;
    prSol[1] =  (prRhs[1] + dt1*(Bx*prRhs[2] + Bz*(-prRhs[0] + prRhs[3]) - By*prRhs[4])*qbym + dtsq*(4*Bx2*prRhs[1] + 4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2] + Bx*(3*Bz*prRhs[4] + By*(prRhs[0] + prRhs[3] - 2*prRhs[5])))*qb2 + dt3*(4*Bx3*prRhs[2] - 2*Bx*(By2 + Bz2)*prRhs[2] + Bz3*(-prRhs[0] + prRhs[3]) - 4*By3*prRhs[4] + 2*By*Bz2*prRhs[4] - By2*Bz*(prRhs[0] - 4*prRhs[3] + 3*prRhs[5]) + Bx2*(2*By*prRhs[4] + Bz*(-4*prRhs[0] + prRhs[3] + 3*prRhs[5])))*qb3 + 2*Bx*By*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[2] =  (prRhs[2] + dt1*(-(Bx*prRhs[1]) + Bz*prRhs[4] + By*(prRhs[0] - prRhs[5]))*qbym + dtsq*(3*By*Bz*prRhs[1] + 4*Bx2*prRhs[2] + By2*prRhs[2] + 4*Bz2*prRhs[2] + Bx*(3*By*prRhs[4] + Bz*(prRhs[0] - 2*prRhs[3] + prRhs[5])))*qb2 + dt3*(-4*Bx3*prRhs[1] + 2*Bx*(By2 + Bz2)*prRhs[1] - 2*By2*Bz*prRhs[4] + 4*Bz3*prRhs[4] + By*Bz2*(prRhs[0] + 3*prRhs[3] - 4*prRhs[5]) + By3*(prRhs[0] - prRhs[5]) - Bx2*(2*Bz*prRhs[4] + By*(-4*prRhs[0] + 3*prRhs[3] + prRhs[5])))*qb3 + 2*Bx*Bz*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[3] =  (prRhs[3] + (-2*Bz*dt1*prRhs[1] + 2*Bx*dt1*prRhs[4])*qbym + dtsq*(2*Bx*By*prRhs[1] + 5*By2*prRhs[3] + Bz2*(2*prRhs[0] + 3*prRhs[3]) + Bz*(-4*Bx*prRhs[2] + 2*By*prRhs[4]) + Bx2*(3*prRhs[3] + 2*prRhs[5]))*qb2 + 2*dt3*(Bx2*(-(Bz*prRhs[1]) + 3*By*prRhs[2]) - Bz*(4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2]) + Bx3*prRhs[4] + Bx*(4*By2*prRhs[4] + Bz2*prRhs[4] + 3*By*Bz*(-prRhs[0] + prRhs[5])))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) + 2*Bx*(2*By2 - Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + 2*By4*prRhs[3] + Bz4*(prRhs[0] + prRhs[3]) + 4*By3*Bz*prRhs[4] - 2*By*Bz3*prRhs[4] + Bx4*(prRhs[3] + prRhs[5]) + By2*Bz2*(prRhs[0] + 3*prRhs[5]) + Bx2*(-2*By*Bz*prRhs[4] + By2*(3*prRhs[0] + prRhs[5]) + Bz2*(prRhs[0] + 2*prRhs[3] + prRhs[5])))*qb4)/d;
    prSol[4] =  (prRhs[4] + dt1*(By*prRhs[1] - Bz*prRhs[2] + Bx*(-prRhs[3] + prRhs[5]))*qbym + dtsq*(3*Bx*Bz*prRhs[1] + Bx2*prRhs[4] + 4*By2*prRhs[4] + 4*Bz2*prRhs[4] + By*(3*Bx*prRhs[2] + Bz*(-2*prRhs[0] + prRhs[3] + prRhs[5])))*qb2 + dt3*(4*By3*prRhs[1] - 2*By*Bz2*prRhs[1] + 2*By2*Bz*prRhs[2] - 4*Bz3*prRhs[2] + Bx2*(-2*By*prRhs[1] + 2*Bz*prRhs[2]) + Bx3*(-prRhs[3] + prRhs[5]) + Bx*(-(Bz2*(3*prRhs[0] + prRhs[3] - 4*prRhs[5])) + By2*(3*prRhs[0] - 4*prRhs[3] + prRhs[5])))*qb3 - 2*By*Bz*dt4*(-6*Bx*(By*prRhs[1] + Bz*prRhs[2]) - 6*By*Bz*prRhs[4] + Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]) + Bx2*(-2*prRhs[0] + prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[5] =  (prRhs[5] + 2*dt1*(By*prRhs[2] - Bx*prRhs[4])*qbym + dtsq*(2*Bx*Bz*prRhs[2] + By*(-4*Bx*prRhs[1] + 2*Bz*prRhs[4]) + 5*Bz2*prRhs[5] + By2*(2*prRhs[0] + 3*prRhs[5]) + Bx2*(2*prRhs[3] + 3*prRhs[5]))*qb2 - 2*dt3*(Bx2*(3*Bz*prRhs[1] - By*prRhs[2]) - By*(3*By*Bz*prRhs[1] + By2*prRhs[2] + 4*Bz2*prRhs[2]) + Bx3*prRhs[4] + Bx*(3*By*Bz*(-prRhs[0] + prRhs[3]) + By2*prRhs[4] + 4*Bz2*prRhs[4]))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 - 2*Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + By2*Bz2*(prRhs[0] + 3*prRhs[3]) - 2*By3*Bz*prRhs[4] + 4*By*Bz3*prRhs[4] + 2*Bz4*prRhs[5] + By4*(prRhs[0] + prRhs[5]) + Bx4*(prRhs[3] + prRhs[5]) + Bx2*(Bz2*(3*prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[3] + 2*prRhs[5])))*qb4)/d;

    for (unsigned i=0; i<6; ++i)
    {
      prTen[6 * n + i] = 2 * prSol[i] - prRhs[i];
    }
  }
  
  return prTen;
}

static std::vector<double>
gkylPressureTensorSrcExactPre(
    MomentSrcData_t *sd, FluidData_t *fd, const double dt, double **ff,
    double *em, double *staticEm)
{
  double Bx = em[BX];
  double By = em[BY];
  double Bz = em[BZ];
  if (sd->hasStatic)
  {
    Bx += staticEm[BX];
    By += staticEm[BY];
    Bz += staticEm[BZ];
  }
  double Bmag = std::sqrt(Bx*Bx + By*By + Bz*Bz);

  unsigned nFluids = sd->nFluids;
  std::vector<double> prTen(6 * nFluids, 0.);
  if (Bmag < 1e-9)  // TODO threshold
    return prTen;

  double bx = Bx / Bmag;
  double by = By / Bmag;
  double bz = Bz / Bmag;

  for (unsigned n = 0; n < nFluids; ++n)
  {
    if (!fd[n].evolve)
      continue;
    double *f = ff[n];
    double Wc = Bmag * fd[n].charge / fd[n].mass;
    double a = -Wc * dt;  // counter-clockwise rotation angle
    double c = std::cos(a);
    double s = std::sin(a);
    double h = 1 - c;  // 2*sin(a/2)^2

    // e1.e1_, e1.e2_, e1.e3_
    double d11 = bx*bx*h + c;
    double d12 = bx*by*h - bz*s;
    double d13 = bx*bz*h + by*s;

    // e2.e1_, e2.e2_, e2.e3_
    double d21 = bx*by*h + bz*s;
    double d22 = by*by*h + c;
    double d23 = by*bz*h - bx*s;
    
    // e3.e1_, e3.e2_, e3.e3_
    double d31 = bx*bz*h - by*s;
    double d32 = by*bz*h + bx*s;  
    double d33 = bz*bz*h + c;

    double pr11 = f[P11] - f[MX]*f[MX]/f[RHO];
    double pr12 = f[P12] - f[MX]*f[MY]/f[RHO];
    double pr13 = f[P13] - f[MX]*f[MZ]/f[RHO];
    double pr22 = f[P22] - f[MY]*f[MY]/f[RHO];
    double pr23 = f[P23] - f[MY]*f[MZ]/f[RHO];
    double pr33 = f[P33] - f[MZ]*f[MZ]/f[RHO];

    prTen[6 *n + 0] = d11*(d11*pr11 + d12*pr12 + d13*pr13)
                    + d12*(d11*pr12 + d12*pr22 + d13*pr23)
                    + d13*(d11*pr13 + d12*pr23 + d13*pr33);
    prTen[6 *n + 1] = d21*(d11*pr11 + d12*pr12 + d13*pr13)
                    + d22*(d11*pr12 + d12*pr22 + d13*pr23)
                    + d23*(d11*pr13 + d12*pr23 + d13*pr33);
    prTen[6 *n + 2] = d31*(d11*pr11 + d12*pr12 + d13*pr13)
                    + d32*(d11*pr12 + d12*pr22 + d13*pr23)
                    + d33*(d11*pr13 + d12*pr23 + d13*pr33);
    prTen[6 *n + 3] = d21*(d21*pr11 + d22*pr12 + d23*pr13)
                    + d22*(d21*pr12 + d22*pr22 + d23*pr23)
                    + d23*(d21*pr13 + d22*pr23 + d23*pr33);
    prTen[6 *n + 4] = d31*(d21*pr11 + d22*pr12 + d23*pr13)
                    + d32*(d21*pr12 + d22*pr22 + d23*pr23)
                    + d33*(d21*pr13 + d22*pr23 + d23*pr33);
    prTen[6 *n + 5] = d31*(d31*pr11 + d32*pr12 + d33*pr13)
                    + d32*(d31*pr12 + d32*pr22 + d33*pr23)
                    + d33*(d31*pr13 + d32*pr23 + d33*pr33);
  }
  
  return prTen;
}
static void
gkylPressureTensorSrcTimeCenteredPost(
    MomentSrcData_t *sd, FluidData_t *fd, const double dt, double **ff,
    const std::vector<double> &prTen)
{
  unsigned nFluids = sd->nFluids;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    if (!fd[n].evolve)
      continue;
    double *f = ff[n];
    f[P11] = f[MX] * f[MX] / f[RHO] + prTen[6 * n + 0];
    f[P12] = f[MX] * f[MY] / f[RHO] + prTen[6 * n + 1];
    f[P13] = f[MX] * f[MZ] / f[RHO] + prTen[6 * n + 2];
    f[P22] = f[MY] * f[MY] / f[RHO] + prTen[6 * n + 3];
    f[P23] = f[MY] * f[MZ] / f[RHO] + prTen[6 * n + 4];
    f[P33] = f[MZ] * f[MZ] / f[RHO] + prTen[6 * n + 5];
  }
}


void
gkylTenMomentSrcTimeCentered(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em,
    double *staticEm)
{
  // update pressure tensor in comoving frame
  std::vector<double> prTen = gkylPressureTensorSrcTimeCenteredPre(
      sd, fd, dt, ff, em, staticEm);

  // update momenta and E field
  gkylMomentSrcTimeCentered(sd, fd, dt, ff, em, staticEm);

  // update pressure tensor in stationary frame
  gkylPressureTensorSrcTimeCenteredPost(sd, fd, dt, ff, prTen);
}


void
gkylTenMomentSrcTimeCenteredDirect(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em,
    double *staticEm)
{
  // update pressure tensor in comoving frame
  std::vector<double> prTen = gkylPressureTensorSrcTimeCenteredDirectPre(
      sd, fd, dt, ff, em, staticEm);

  // update momenta and E field
  gkylMomentSrcTimeCenteredDirect(sd, fd, dt, ff, em, staticEm);

  // update pressure tensor in stationary frame
  gkylPressureTensorSrcTimeCenteredPost(sd, fd, dt, ff, prTen);
}


void
gkylTenMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  // update pressure tensor in comoving frame
  std::vector<double> prTen = gkylPressureTensorSrcTimeCenteredDirectPre(
      sd, fd, dt, ff, em, staticEm);

  // update momenta and E field
  gkylMomentSrcTimeCenteredDirect2(sd, fd, dt, ff, em, staticEm);

  // update pressure tensor in stationary frame
  gkylPressureTensorSrcTimeCenteredPost(sd, fd, dt, ff, prTen);
}


void
gkylTenMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  // update pressure tensor in comoving frame
  std::vector<double> prTen = gkylPressureTensorSrcExactPre(
      sd, fd, dt, ff, em, staticEm);

  // update momenta and E field
  gkylMomentSrcExact(sd, fd, dt, ff, em, staticEm);

  // update pressure tensor in stationary frame
  gkylPressureTensorSrcTimeCenteredPost(sd, fd, dt, ff, prTen);
}

