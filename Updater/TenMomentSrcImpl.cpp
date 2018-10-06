// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
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
static const unsigned ER = 4;

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

static double sq(double x) { return x*x; }

static const int COL_PIV_HOUSEHOLDER_QR = 0;
static const int PARTIAL_PIV_LU = 1;

#define fidx(n, c) (3 * (n) + (c))
#define eidx(c) (3 * nFluids + (c))

void
gkylTenMomentSrcRk3(TenMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em)
{
}

void
gkylTenMomentSrcTimeCentered(TenMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  unsigned nFluids = sd->nFluids;
  double dt1 = 0.5 * dt;
  double dt2 = 0.5 * dt / sd->epsilon0;
  
  std::vector<double> zeros(6, 0);
  if (!(sd->hasStatic))
  {
    staticEm = zeros.data();
  }

  // for momentum and electric field update
  Eigen::MatrixXd lhs;
  lhs = Eigen::MatrixXd::Constant(3*nFluids+3, 3*nFluids+3, 0.0);
  Eigen::VectorXd rhs(3*nFluids+3);

  // for pressure update
  Eigen::MatrixXd prLhs;
  prLhs = Eigen::MatrixXd::Constant(6, 6, 0.0);
  Eigen::VectorXd prRhs(6);
  // updated pressure tensor
  std::vector<double> prTen(6*nFluids);

  double bx = em[BX]+staticEm[BX];
  double by = em[BY]+staticEm[BY];
  double bz = em[BZ]+staticEm[BZ];

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


  // compute increments from pressure tensor terms
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;

    // assemble LHS terms
    prLhs(0,0) = 1;
    prLhs(0,1) = -2*dt1*qbym*bz;
    prLhs(0,2) = 2*dt1*qbym*by;
    prLhs(0,3) = 0;
    prLhs(0,4) = 0;
    prLhs(0,5) = 0;
    prLhs(1,0) = dt1*qbym*bz;
    prLhs(1,1) = 1;
    prLhs(1,2) = -dt1*qbym*bx;
    prLhs(1,3) = -dt1*qbym*bz;
    prLhs(1,4) = dt1*qbym*by;
    prLhs(1,5) = 0;
    prLhs(2,0) = -dt1*qbym*by;
    prLhs(2,1) = dt1*qbym*bx;
    prLhs(2,2) = 1;
    prLhs(2,3) = 0;
    prLhs(2,4) = -dt1*qbym*bz;
    prLhs(2,5) = dt1*qbym*by;
    prLhs(3,0) = 0;
    prLhs(3,1) = 2*dt1*qbym*bz;
    prLhs(3,2) = 0;
    prLhs(3,3) = 1;
    prLhs(3,4) = -2*dt1*qbym*bx;
    prLhs(3,5) = 0;
    prLhs(4,0) = 0;
    prLhs(4,1) = -dt1*qbym*by;
    prLhs(4,2) = dt1*qbym*bz;
    prLhs(4,3) = dt1*qbym*bx;
    prLhs(4,4) = 1;
    prLhs(4,5) = -dt1*qbym*bx;
    prLhs(5,0) = 0;
    prLhs(5,1) = 0;
    prLhs(5,2) = -2*dt1*qbym*by;
    prLhs(5,3) = 0;
    prLhs(5,4) = 2*dt1*qbym*bx;
    prLhs(5,5) = 1;

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
      prTen[6*n+i] = 2*prSol[i]-prRhs[i];
  }

  //------------> update solution for fluids
  double chargeDens = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;

    chargeDens += qbym*f[RHO];

    // update momentum
    f[MX] = 2*sol(fidx(n,X))/qbym - f[MX];
    f[MY] = 2*sol(fidx(n,Y))/qbym - f[MY];
    f[MZ] = 2*sol(fidx(n,Z))/qbym - f[MZ];
    
    // update pressure tensor
    f[P11] = f[MX]*f[MX]/f[RHO] + prTen[6*n+0];
    f[P12] = f[MX]*f[MY]/f[RHO] + prTen[6*n+1];
    f[P13] = f[MX]*f[MZ]/f[RHO] + prTen[6*n+2];
    f[P22] = f[MY]*f[MY]/f[RHO] + prTen[6*n+3];
    f[P23] = f[MY]*f[MZ]/f[RHO] + prTen[6*n+4];
    f[P33] = f[MZ]*f[MZ]/f[RHO] + prTen[6*n+5];
  }

  //------------> update electric field
  em[EX] = 2*sol(eidx(X)) - em[EX];
  em[EY] = 2*sol(eidx(Y)) - em[EY];
  em[EZ] = 2*sol(eidx(Z)) - em[EZ];

  //------------> update correction potential
  double crhoc = sd->chi_e*chargeDens/sd->epsilon0;
  em[PHIE] += dt*crhoc;
}
