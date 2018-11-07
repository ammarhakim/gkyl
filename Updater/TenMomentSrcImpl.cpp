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

  double Bx = em[BX]+staticEm[BX];
  double By = em[BY]+staticEm[BY];
  double Bz = em[BZ]+staticEm[BZ];

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
    prLhs(0,1) = -2*dt1*qbym*Bz;
    prLhs(0,2) = 2*dt1*qbym*By;
    prLhs(0,3) = 0;
    prLhs(0,4) = 0;
    prLhs(0,5) = 0;
    prLhs(1,0) = dt1*qbym*Bz;
    prLhs(1,1) = 1;
    prLhs(1,2) = -dt1*qbym*Bx;
    prLhs(1,3) = -dt1*qbym*Bz;
    prLhs(1,4) = dt1*qbym*By;
    prLhs(1,5) = 0;
    prLhs(2,0) = -dt1*qbym*By;
    prLhs(2,1) = dt1*qbym*Bx;
    prLhs(2,2) = 1;
    prLhs(2,3) = 0;
    prLhs(2,4) = -dt1*qbym*Bz;
    prLhs(2,5) = dt1*qbym*By;
    prLhs(3,0) = 0;
    prLhs(3,1) = 2*dt1*qbym*Bz;
    prLhs(3,2) = 0;
    prLhs(3,3) = 1;
    prLhs(3,4) = -2*dt1*qbym*Bx;
    prLhs(3,5) = 0;
    prLhs(4,0) = 0;
    prLhs(4,1) = -dt1*qbym*By;
    prLhs(4,2) = dt1*qbym*Bz;
    prLhs(4,3) = dt1*qbym*Bx;
    prLhs(4,4) = 1;
    prLhs(4,5) = -dt1*qbym*Bx;
    prLhs(5,0) = 0;
    prLhs(5,1) = 0;
    prLhs(5,2) = -2*dt1*qbym*By;
    prLhs(5,3) = 0;
    prLhs(5,4) = 2*dt1*qbym*Bx;
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
    {
      prTen[6*n+i] = 2*prSol[i]-prRhs[i];
    }
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

void
gkylTenMomentSrcAnalytic(TenMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  // based on Smithe (2007) with corrections
  unsigned nFluids = sd->nFluids;
  double epsilon0 = sd->epsilon0;
  
  std::vector<double> zeros(6, 0);
  if (!(sd->hasStatic))
  {
    staticEm = zeros.data();
  }

  std::vector<double> lambda(nFluids);
  std::vector<double> omega(nFluids);
  //      std::vector<double> K(3); 
  //      double wp2 = 0.0;
  double gamma2 = 0.0;
  double delta = 0.0;
  double theta = 0.0; 
  double Bx = (em[BX] + staticEm[BX]);
  double By = (em[BY] + staticEm[BY]);
  double Bz = (em[BZ] + staticEm[BZ]);
  Eigen::Vector3d B; 
  Eigen::Vector3d E;
  Eigen::Vector3d F(0,0,0); 
  Eigen::Vector3d K(0,0,0);// the k vector used to update the implicit solution in Smithe(2007)
  B(0) = Bx; B(1) = By; B(2) = Bz;
  E(0) = em[EX]; E(1) = em[EY]; E(2) = em[EZ];
  double babs = std::sqrt(Bx*Bx + By*By + Bz*Bz);

  //------------> fill elements corresponding to each fluid
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;
    double qbym2 = sq(qbym);

    lambda[n] = 1.00; //collisional parameter
    omega[n] = qbym*dt; // omega/B
    double wp2 = f[RHO]*qbym2/epsilon0*dt*dt;
    //          wp2 += f[RHO]*qbym2/epsilon0*dt*dt;
    gamma2 += wp2*qbym2*dt*dt/(1+0.25*qbym2*babs*babs*dt*dt);
    delta += wp2*qbym*dt/(1+0.25*qbym2*babs*babs*dt*dt);
    theta += wp2/(1+qbym2*dt*dt/4*babs*babs);
    Eigen::Vector3d j(f[MX]*qbym, f[MY]*qbym, f[MZ]*qbym);
    double bdotj = B.dot(j);
    Eigen::Vector3d bcrossj = B.cross(j);

    K -= dt/2.0 *(1.0 + lambda[n])*(1.0*j + 0.25*omega[n]*omega[n]*B * bdotj - 0.5*omega[n]*bcrossj )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
    F -= dt/2.0*(1/4.0*omega[n]*omega[n]*f[RHO]*qbym2/epsilon0*dt/2.0*B*B.dot(E) + 1/2.0 * f[RHO]*qbym2/epsilon0*dt*E - 1/2.0*omega[n]*f[RHO]*qbym2/epsilon0*dt*B.cross(E)/2.0 )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
  }

  F = F + E;

  double denom = ((1+theta/4)*(1+theta/4) + delta*delta*babs*babs/64.0) * (1 + theta/4 + gamma2/16*babs*babs);        
  double coeff_e = ((1+theta/4)*(1+theta/4) + gamma2*babs*babs/16*(1+theta/4));
  //        double coeff_k = 1.0/denom*0; 
  double coeff_dot = (1.0/64.0*delta*delta - 1.0/16.0*gamma2*(1+theta/4));
  double coeff_cross = 0.125*delta*(1 + gamma2*babs*babs/16 + theta/4);
  Eigen::Vector3d fk = 1.0*F + K/epsilon0; // vector holding sum of K + 2 epsilon0 E used in many calculations
  //        fk(0) = K[0] + em[EX]*2*epsilon0;
  //        fk(1) = K[1] + em[EY]*2*epsilon0; 
  //        fk(2) = K[2] + em[EZ]*2*epsilon0; 
  em[EX] = (coeff_e*fk[EX] + coeff_dot*B(0)*B.dot(fk) + coeff_cross*(B.cross(fk)(0)) )/denom; 
  em[EY] = (coeff_e*fk[EY] + coeff_dot*B(1)*B.dot(fk) + coeff_cross*(B.cross(fk)(1)) )/denom; 
  em[EZ] = (coeff_e*fk[EZ] + coeff_dot*B(2)*B.dot(fk) + coeff_cross*(B.cross(fk)(2)) )/denom; 

  // update the stored E field to E_n+1/2
  E(0) = (em[EX] + E(0))/2.0;
  E(1) = (em[EY] + E(1))/2.0;
  E(2) = (em[EZ] + E(2))/2.0;

  // pressure tensor update
  std::vector<double> prTen(6*nFluids);
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
  std::vector<double> qbym(nFluids); 
  for (unsigned n = 0; n < nFluids; ++n)
    qbym[n] = fd[n].charge/fd[n].mass;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];

    prRhs[0] = f[P11] - f[MX] * f[MX] / f[RHO];
    prRhs[1] = f[P12] - f[MX] * f[MY] / f[RHO];
    prRhs[2] = f[P13] - f[MX] * f[MZ] / f[RHO];
    prRhs[3] = f[P22] - f[MY] * f[MY] / f[RHO];
    prRhs[4] = f[P23] - f[MY] * f[MZ] / f[RHO];
    prRhs[5] = f[P33] - f[MZ] * f[MZ] / f[RHO];

    double qb2 = qbym[n] * qbym[n];
    double qb3 = qbym[n] * qb2;
    double qb4 = qbym[n] * qb3;
    double d = 1 + 5*(Bx2 + By2 + Bz2)*dtsq*qb2 + 4*(Bx2 + By2 + Bz2)*(Bx2 + By2 + Bz2)*dt4*qb4;

    prSol[0] = (prRhs[0] + 2*dt1*(Bz*prRhs[1] - By*prRhs[2])*qbym[n] + dtsq*(5*Bx2*prRhs[0] + 2*Bx*(By*prRhs[1] + Bz*prRhs[2]) + Bz2*(3*prRhs[0] + 2*prRhs[3]) - 4*By*Bz*prRhs[4] + By2*(3*prRhs[0] + 2*prRhs[5]))*qb2 + 2*dt3*(4*Bx2*(Bz*prRhs[1] - By*prRhs[2]) - (By2 + Bz2)*(-(Bz*prRhs[1]) + By*prRhs[2]) - 3*Bx*(By2*prRhs[4] - Bz2*prRhs[4] + By*Bz*(-prRhs[3] + prRhs[5])))*qb3 + 2*dt4*(2*Bx4*prRhs[0] + 4*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 + Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + (By2 + Bz2)*(Bz2*(prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[5])) + Bx2*(4*By*Bz*prRhs[4] + By2*(3*prRhs[3] + prRhs[5]) + Bz2*(prRhs[3] + 3*prRhs[5])))*qb4)/d;
    prSol[1] =  (prRhs[1] + dt1*(Bx*prRhs[2] + Bz*(-prRhs[0] + prRhs[3]) - By*prRhs[4])*qbym[n] + dtsq*(4*Bx2*prRhs[1] + 4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2] + Bx*(3*Bz*prRhs[4] + By*(prRhs[0] + prRhs[3] - 2*prRhs[5])))*qb2 + dt3*(4*Bx3*prRhs[2] - 2*Bx*(By2 + Bz2)*prRhs[2] + Bz3*(-prRhs[0] + prRhs[3]) - 4*By3*prRhs[4] + 2*By*Bz2*prRhs[4] - By2*Bz*(prRhs[0] - 4*prRhs[3] + 3*prRhs[5]) + Bx2*(2*By*prRhs[4] + Bz*(-4*prRhs[0] + prRhs[3] + 3*prRhs[5])))*qb3 + 2*Bx*By*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[2] =  (prRhs[2] + dt1*(-(Bx*prRhs[1]) + Bz*prRhs[4] + By*(prRhs[0] - prRhs[5]))*qbym[n] + dtsq*(3*By*Bz*prRhs[1] + 4*Bx2*prRhs[2] + By2*prRhs[2] + 4*Bz2*prRhs[2] + Bx*(3*By*prRhs[4] + Bz*(prRhs[0] - 2*prRhs[3] + prRhs[5])))*qb2 + dt3*(-4*Bx3*prRhs[1] + 2*Bx*(By2 + Bz2)*prRhs[1] - 2*By2*Bz*prRhs[4] + 4*Bz3*prRhs[4] + By*Bz2*(prRhs[0] + 3*prRhs[3] - 4*prRhs[5]) + By3*(prRhs[0] - prRhs[5]) - Bx2*(2*Bz*prRhs[4] + By*(-4*prRhs[0] + 3*prRhs[3] + prRhs[5])))*qb3 + 2*Bx*Bz*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[3] =  (prRhs[3] + (-2*Bz*dt1*prRhs[1] + 2*Bx*dt1*prRhs[4])*qbym[n] + dtsq*(2*Bx*By*prRhs[1] + 5*By2*prRhs[3] + Bz2*(2*prRhs[0] + 3*prRhs[3]) + Bz*(-4*Bx*prRhs[2] + 2*By*prRhs[4]) + Bx2*(3*prRhs[3] + 2*prRhs[5]))*qb2 + 2*dt3*(Bx2*(-(Bz*prRhs[1]) + 3*By*prRhs[2]) - Bz*(4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2]) + Bx3*prRhs[4] + Bx*(4*By2*prRhs[4] + Bz2*prRhs[4] + 3*By*Bz*(-prRhs[0] + prRhs[5])))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) + 2*Bx*(2*By2 - Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + 2*By4*prRhs[3] + Bz4*(prRhs[0] + prRhs[3]) + 4*By3*Bz*prRhs[4] - 2*By*Bz3*prRhs[4] + Bx4*(prRhs[3] + prRhs[5]) + By2*Bz2*(prRhs[0] + 3*prRhs[5]) + Bx2*(-2*By*Bz*prRhs[4] + By2*(3*prRhs[0] + prRhs[5]) + Bz2*(prRhs[0] + 2*prRhs[3] + prRhs[5])))*qb4)/d;
    prSol[4] =  (prRhs[4] + dt1*(By*prRhs[1] - Bz*prRhs[2] + Bx*(-prRhs[3] + prRhs[5]))*qbym[n] + dtsq*(3*Bx*Bz*prRhs[1] + Bx2*prRhs[4] + 4*By2*prRhs[4] + 4*Bz2*prRhs[4] + By*(3*Bx*prRhs[2] + Bz*(-2*prRhs[0] + prRhs[3] + prRhs[5])))*qb2 + dt3*(4*By3*prRhs[1] - 2*By*Bz2*prRhs[1] + 2*By2*Bz*prRhs[2] - 4*Bz3*prRhs[2] + Bx2*(-2*By*prRhs[1] + 2*Bz*prRhs[2]) + Bx3*(-prRhs[3] + prRhs[5]) + Bx*(-(Bz2*(3*prRhs[0] + prRhs[3] - 4*prRhs[5])) + By2*(3*prRhs[0] - 4*prRhs[3] + prRhs[5])))*qb3 - 2*By*Bz*dt4*(-6*Bx*(By*prRhs[1] + Bz*prRhs[2]) - 6*By*Bz*prRhs[4] + Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]) + Bx2*(-2*prRhs[0] + prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[5] =  (prRhs[5] + 2*dt1*(By*prRhs[2] - Bx*prRhs[4])*qbym[n] + dtsq*(2*Bx*Bz*prRhs[2] + By*(-4*Bx*prRhs[1] + 2*Bz*prRhs[4]) + 5*Bz2*prRhs[5] + By2*(2*prRhs[0] + 3*prRhs[5]) + Bx2*(2*prRhs[3] + 3*prRhs[5]))*qb2 - 2*dt3*(Bx2*(3*Bz*prRhs[1] - By*prRhs[2]) - By*(3*By*Bz*prRhs[1] + By2*prRhs[2] + 4*Bz2*prRhs[2]) + Bx3*prRhs[4] + Bx*(3*By*Bz*(-prRhs[0] + prRhs[3]) + By2*prRhs[4] + 4*Bz2*prRhs[4]))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 - 2*Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + By2*Bz2*(prRhs[0] + 3*prRhs[3]) - 2*By3*Bz*prRhs[4] + 4*By*Bz3*prRhs[4] + 2*Bz4*prRhs[5] + By4*(prRhs[0] + prRhs[5]) + Bx4*(prRhs[3] + prRhs[5]) + Bx2*(Bz2*(3*prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[3] + 2*prRhs[5])))*qb4)/d;

    for (unsigned i=0; i<6; ++i)
    {
      prTen[6*n+i] = 2*prSol[i]-prRhs[i];
    }
  }

  //------------> update solution for fluids
  double chargeDens = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;
    double qbym2 = sq(qbym);

    chargeDens += qbym*f[RHO];

    double j_coeff = (lambda[n] - 0.25*omega[n]*omega[n]*babs*babs)/(1.0+0.25*omega[n]*omega[n]*babs*babs);
    double f_coeff = f[RHO]*qbym2/epsilon0*dt/(1.0+0.25*omega[n]*omega[n]*babs*babs);
    
    double dot_coeff = 0.25*omega[n]*omega[n]/(1.0+0.25*omega[n]*omega[n]*babs*babs);
    double cross_coeff = omega[n]/2.0/(1.0+0.25*omega[n]*omega[n]*babs*babs);
    // update E with the new solution
    Eigen::Vector3d j(f[MX]*qbym,f[MY]*qbym,f[MZ]*qbym);
    Eigen::Vector3d jf = (1.0+lambda[n])*j + E*epsilon0*dt*f[RHO]*qbym2/epsilon0; // There is an error in smithe's paper here: F in equation (15) must be multiplied by a factor of wp^2 dt or the units are wrong.
    Eigen::Vector3d solj = j_coeff*j + f_coeff*E*epsilon0 + dot_coeff*B*B.dot(jf) - cross_coeff*B.cross(jf); 

    f[MX] = solj(0)/qbym;
    f[MY] = solj(1)/qbym;
    f[MZ] = solj(2)/qbym;

    f[P11] = f[MX] * f[MX] / f[RHO] + prTen[6 * n + 0];
    f[P12] = f[MX] * f[MY] / f[RHO] + prTen[6 * n + 1];
    f[P13] = f[MX] * f[MZ] / f[RHO] + prTen[6 * n + 2];
    f[P22] = f[MY] * f[MY] / f[RHO] + prTen[6 * n + 3];
    f[P23] = f[MY] * f[MZ] / f[RHO] + prTen[6 * n + 4];
    f[P33] = f[MZ] * f[MZ] / f[RHO] + prTen[6 * n + 5];
  }


  //------------> update correction potential
  double crhoc = sd->chi_e*chargeDens/sd->epsilon0;
  em[PHIE] += dt*crhoc;
}


void
gkylTenMomentSrcAnalytic2(TenMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  // based on Smithe (2007) with corrections
  // but using Hakim (2019) notations
  unsigned nFluids = sd->nFluids;
  double epsilon0 = sd->epsilon0;

  std::vector<double> zeros(6, 0);
  if (!(sd->hasStatic))
  {
    staticEm = zeros.data();
  }

  double Bx = (em[BX] + staticEm[BX]);
  double By = (em[BY] + staticEm[BY]);
  double Bz = (em[BZ] + staticEm[BZ]);
  double Bmag = std::sqrt(Bx*Bx + By*By + Bz*Bz);
  Eigen::Vector3d b(Bx / Bmag, By / Bmag, Bz / Bmag);

  std::vector<double> qbym(nFluids); 
  std::vector<Eigen::Vector3d> J(nFluids); 
  double w02 = 0.;
  double gam2 = 0.;
  double delta = 0.;
  Eigen::Vector3d K(0., 0., 0.);
  std::vector<double> Wc_dt(nFluids); 
  std::vector<double> wp_dt2(nFluids); 
  for (unsigned n=0; n < nFluids; ++n)
  {
    qbym[n] = fd[n].charge / fd[n].mass;
    double *f = ff[n];
    J[n][0] = f[MX] * qbym[n];
    J[n][1] = f[MY] * qbym[n];
    J[n][2] = f[MZ] * qbym[n];
    if (!fd[n].evolve)
      continue;
    Wc_dt[n] = qbym[n] * Bmag * dt;
    wp_dt2[n] = f[RHO] * sq(qbym[n]) / epsilon0 * sq(dt);
    double tmp = 1. + sq(Wc_dt[n]) / 4.;
    w02 += wp_dt2[n] / tmp;
    gam2 += wp_dt2[n] * sq(Wc_dt[n]) / tmp;
    delta += wp_dt2[n] * Wc_dt[n] / tmp;
    K -= dt / tmp * (J[n] + sq(Wc_dt[n] / 2.) * b * b.dot(J[n]) - (Wc_dt[n] / 2.) * b.cross(J[n]));
  }
  double Delta2 = sq(delta) / (1. + w02 / 4.);

  Eigen::Vector3d F(em[EX] * epsilon0, em[EY] * epsilon0, em[EZ] * epsilon0);
  Eigen::Vector3d F_halfK = F + 0.5 * K;
  for (unsigned n=0; n < nFluids; ++n)
  {
    if (fd[n].evolve)
      continue;
    F_halfK -= (0.5 * dt) * J[n];
  }

  double tmp = 1. / (1. + w02 / 4. + Delta2 / 64.);
  Eigen::Vector3d Fbar = tmp * (
      F_halfK 
      + ((Delta2 / 64. - gam2 / 16.) / (1. + w02 / 4. + gam2 / 16.)) * b * b.dot(F_halfK)
      + (delta / 8. / (1. + w02 / 4.)) * b.cross(F_halfK)
      );

  Eigen::Vector3d F_new = 2. * Fbar - F;
  em[EX] = F_new[0] / epsilon0;
  em[EY] = F_new[1] / epsilon0;
  em[EZ] = F_new[2] / epsilon0;

  // pressure tensor update
  std::vector<double> prTen(6*nFluids);
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
    double *f = ff[n];

    prRhs[0] = f[P11] - f[MX] * f[MX] / f[RHO];
    prRhs[1] = f[P12] - f[MX] * f[MY] / f[RHO];
    prRhs[2] = f[P13] - f[MX] * f[MZ] / f[RHO];
    prRhs[3] = f[P22] - f[MY] * f[MY] / f[RHO];
    prRhs[4] = f[P23] - f[MY] * f[MZ] / f[RHO];
    prRhs[5] = f[P33] - f[MZ] * f[MZ] / f[RHO];

    double qb2 = qbym[n] * qbym[n];
    double qb3 = qbym[n] * qb2;
    double qb4 = qbym[n] * qb3;
    double d = 1 + 5*(Bx2 + By2 + Bz2)*dtsq*qb2 + 4*(Bx2 + By2 + Bz2)*(Bx2 + By2 + Bz2)*dt4*qb4;

    prSol[0] = (prRhs[0] + 2*dt1*(Bz*prRhs[1] - By*prRhs[2])*qbym[n] + dtsq*(5*Bx2*prRhs[0] + 2*Bx*(By*prRhs[1] + Bz*prRhs[2]) + Bz2*(3*prRhs[0] + 2*prRhs[3]) - 4*By*Bz*prRhs[4] + By2*(3*prRhs[0] + 2*prRhs[5]))*qb2 + 2*dt3*(4*Bx2*(Bz*prRhs[1] - By*prRhs[2]) - (By2 + Bz2)*(-(Bz*prRhs[1]) + By*prRhs[2]) - 3*Bx*(By2*prRhs[4] - Bz2*prRhs[4] + By*Bz*(-prRhs[3] + prRhs[5])))*qb3 + 2*dt4*(2*Bx4*prRhs[0] + 4*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 + Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + (By2 + Bz2)*(Bz2*(prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[5])) + Bx2*(4*By*Bz*prRhs[4] + By2*(3*prRhs[3] + prRhs[5]) + Bz2*(prRhs[3] + 3*prRhs[5])))*qb4)/d;
    prSol[1] =  (prRhs[1] + dt1*(Bx*prRhs[2] + Bz*(-prRhs[0] + prRhs[3]) - By*prRhs[4])*qbym[n] + dtsq*(4*Bx2*prRhs[1] + 4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2] + Bx*(3*Bz*prRhs[4] + By*(prRhs[0] + prRhs[3] - 2*prRhs[5])))*qb2 + dt3*(4*Bx3*prRhs[2] - 2*Bx*(By2 + Bz2)*prRhs[2] + Bz3*(-prRhs[0] + prRhs[3]) - 4*By3*prRhs[4] + 2*By*Bz2*prRhs[4] - By2*Bz*(prRhs[0] - 4*prRhs[3] + 3*prRhs[5]) + Bx2*(2*By*prRhs[4] + Bz*(-4*prRhs[0] + prRhs[3] + 3*prRhs[5])))*qb3 + 2*Bx*By*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[2] =  (prRhs[2] + dt1*(-(Bx*prRhs[1]) + Bz*prRhs[4] + By*(prRhs[0] - prRhs[5]))*qbym[n] + dtsq*(3*By*Bz*prRhs[1] + 4*Bx2*prRhs[2] + By2*prRhs[2] + 4*Bz2*prRhs[2] + Bx*(3*By*prRhs[4] + Bz*(prRhs[0] - 2*prRhs[3] + prRhs[5])))*qb2 + dt3*(-4*Bx3*prRhs[1] + 2*Bx*(By2 + Bz2)*prRhs[1] - 2*By2*Bz*prRhs[4] + 4*Bz3*prRhs[4] + By*Bz2*(prRhs[0] + 3*prRhs[3] - 4*prRhs[5]) + By3*(prRhs[0] - prRhs[5]) - Bx2*(2*Bz*prRhs[4] + By*(-4*prRhs[0] + 3*prRhs[3] + prRhs[5])))*qb3 + 2*Bx*Bz*dt4*(6*Bx*(By*prRhs[1] + Bz*prRhs[2]) + 6*By*Bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[3] =  (prRhs[3] + (-2*Bz*dt1*prRhs[1] + 2*Bx*dt1*prRhs[4])*qbym[n] + dtsq*(2*Bx*By*prRhs[1] + 5*By2*prRhs[3] + Bz2*(2*prRhs[0] + 3*prRhs[3]) + Bz*(-4*Bx*prRhs[2] + 2*By*prRhs[4]) + Bx2*(3*prRhs[3] + 2*prRhs[5]))*qb2 + 2*dt3*(Bx2*(-(Bz*prRhs[1]) + 3*By*prRhs[2]) - Bz*(4*By2*prRhs[1] + Bz2*prRhs[1] + 3*By*Bz*prRhs[2]) + Bx3*prRhs[4] + Bx*(4*By2*prRhs[4] + Bz2*prRhs[4] + 3*By*Bz*(-prRhs[0] + prRhs[5])))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) + 2*Bx*(2*By2 - Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + 2*By4*prRhs[3] + Bz4*(prRhs[0] + prRhs[3]) + 4*By3*Bz*prRhs[4] - 2*By*Bz3*prRhs[4] + Bx4*(prRhs[3] + prRhs[5]) + By2*Bz2*(prRhs[0] + 3*prRhs[5]) + Bx2*(-2*By*Bz*prRhs[4] + By2*(3*prRhs[0] + prRhs[5]) + Bz2*(prRhs[0] + 2*prRhs[3] + prRhs[5])))*qb4)/d;
    prSol[4] =  (prRhs[4] + dt1*(By*prRhs[1] - Bz*prRhs[2] + Bx*(-prRhs[3] + prRhs[5]))*qbym[n] + dtsq*(3*Bx*Bz*prRhs[1] + Bx2*prRhs[4] + 4*By2*prRhs[4] + 4*Bz2*prRhs[4] + By*(3*Bx*prRhs[2] + Bz*(-2*prRhs[0] + prRhs[3] + prRhs[5])))*qb2 + dt3*(4*By3*prRhs[1] - 2*By*Bz2*prRhs[1] + 2*By2*Bz*prRhs[2] - 4*Bz3*prRhs[2] + Bx2*(-2*By*prRhs[1] + 2*Bz*prRhs[2]) + Bx3*(-prRhs[3] + prRhs[5]) + Bx*(-(Bz2*(3*prRhs[0] + prRhs[3] - 4*prRhs[5])) + By2*(3*prRhs[0] - 4*prRhs[3] + prRhs[5])))*qb3 - 2*By*Bz*dt4*(-6*Bx*(By*prRhs[1] + Bz*prRhs[2]) - 6*By*Bz*prRhs[4] + Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]) + Bx2*(-2*prRhs[0] + prRhs[3] + prRhs[5]))*qb4)/d;
    prSol[5] =  (prRhs[5] + 2*dt1*(By*prRhs[2] - Bx*prRhs[4])*qbym[n] + dtsq*(2*Bx*Bz*prRhs[2] + By*(-4*Bx*prRhs[1] + 2*Bz*prRhs[4]) + 5*Bz2*prRhs[5] + By2*(2*prRhs[0] + 3*prRhs[5]) + Bx2*(2*prRhs[3] + 3*prRhs[5]))*qb2 - 2*dt3*(Bx2*(3*Bz*prRhs[1] - By*prRhs[2]) - By*(3*By*Bz*prRhs[1] + By2*prRhs[2] + 4*Bz2*prRhs[2]) + Bx3*prRhs[4] + Bx*(3*By*Bz*(-prRhs[0] + prRhs[3]) + By2*prRhs[4] + 4*Bz2*prRhs[4]))*qb3 + 2*dt4*(-2*Bx3*(By*prRhs[1] + Bz*prRhs[2]) - 2*Bx*(By2 - 2*Bz2)*(By*prRhs[1] + Bz*prRhs[2]) + By2*Bz2*(prRhs[0] + 3*prRhs[3]) - 2*By3*Bz*prRhs[4] + 4*By*Bz3*prRhs[4] + 2*Bz4*prRhs[5] + By4*(prRhs[0] + prRhs[5]) + Bx4*(prRhs[3] + prRhs[5]) + Bx2*(Bz2*(3*prRhs[0] + prRhs[3]) - 2*By*Bz*prRhs[4] + By2*(prRhs[0] + prRhs[3] + 2*prRhs[5])))*qb4)/d;

    for (unsigned i=0; i<6; ++i)
    {
      prTen[6 * n + i] = 2 * prSol[i] - prRhs[i];
    }
  }

  double chargeDens = 0.0;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    chargeDens += qbym[n] * f[RHO];
    if (!fd[n].evolve)
      continue;
    Eigen::Vector3d Jstar = J[n] + Fbar * (wp_dt2[n] / dt / 2.);
    Eigen::Vector3d J_new = 2. * (Jstar + sq(Wc_dt[n] / 2.) * b * b.dot(Jstar) - (Wc_dt[n] / 2.) * b.cross(Jstar)) / (1. + sq(Wc_dt[n] / 2.)) - J[n];

    f[MX] = J_new[0] / qbym[n];
    f[MY] = J_new[1] / qbym[n];
    f[MZ] = J_new[2] / qbym[n];

    f[P11] = f[MX] * f[MX] / f[RHO] + prTen[6 * n + 0];
    f[P12] = f[MX] * f[MY] / f[RHO] + prTen[6 * n + 1];
    f[P13] = f[MX] * f[MZ] / f[RHO] + prTen[6 * n + 2];
    f[P22] = f[MY] * f[MY] / f[RHO] + prTen[6 * n + 3];
    f[P23] = f[MY] * f[MZ] / f[RHO] + prTen[6 * n + 4];
    f[P33] = f[MZ] * f[MZ] / f[RHO] + prTen[6 * n + 5];
  }

  //------------> update correction potential
  double crhoc = sd->chi_e * chargeDens/sd->epsilon0;
  em[PHIE] += dt * crhoc;
}
