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
  std::vector<double> keOld(nFluids);
  if (sd->hasPressure)
  {
    for (unsigned n=0; n<nFluids; ++n)
      keOld[n] = 0.5*(sq(ff[n][MX]) + sq(ff[n][MY]) + sq(ff[n][MZ]))/ff[n][RHO];
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
      ff[n][ER] += keNew-keOld[n];
    }
  }  
}

void
gkylFiveMomentSrcTimeCentered(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
{
  unsigned nFluids = sd->nFluids;
  double dt1 = 0.5 * dt;
  double dt2 = 0.5 * dt / sd->epsilon0;
  
  std::vector<double> zeros(6, 0);
  if (!(sd->hasStatic))
  {
    staticEm = zeros.data();
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

void
gkylFiveMomentSrcAnalytic(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
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
  double bx = (em[BX] + staticEm[BX]);
  double by = (em[BY] + staticEm[BY]);
  double bz = (em[BZ] + staticEm[BZ]);
  Eigen::Vector3d B; 
  Eigen::Vector3d E;
  Eigen::Vector3d F(0,0,0); 
  Eigen::Vector3d K(0,0,0);// the k vector used to update the implicit solution in Smithe(2007)
  B(0) = bx; B(1) = by; B(2) = bz;
  E(0) = em[EX]; E(1) = em[EY]; E(2) = em[EZ];
  double babs = std::sqrt(bx*bx + by*by + bz*bz);

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

  //------------> update solution for fluids
  double chargeDens = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  {
    double *f = ff[n];
    double qbym = fd[n].charge/fd[n].mass;
    double qbym2 = sq(qbym);

    chargeDens += qbym*f[RHO];

    double keOld = 0.5*(sq(f[MX]) + sq(f[MY]) + sq(f[MZ]))/f[RHO];

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

    if (sd->hasPressure)
    {
      double keNew = 0.5*(sq(f[MX]) + sq(f[MY]) + sq(f[MZ]))/f[RHO];
      f[ER] += keNew-keOld;
    }
  }

  //------------> update correction potential
  double crhoc = sd->chi_e*chargeDens/sd->epsilon0;
  em[PHIE] += dt*crhoc;
}

void
gkylFiveMomentSrcAnalytic2(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
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

  double w02 = 0.;
  double gam2 = 0.;
  double delta = 0.;
  Eigen::Vector3d K(0., 0., 0.);
  std::vector<double> qbym(nFluids); 
  std::vector<double> Wc_dt(nFluids); 
  std::vector<double> wp_dt2(nFluids); 
  for (unsigned n=0; n < nFluids; ++n)
  {
    double *f = ff[n];
    qbym[n] = fd[n].charge / fd[n].mass;
    Wc_dt[n] = qbym[n] * Bmag * dt;
    wp_dt2[n] = f[RHO] * sq(qbym[n]) / epsilon0 * sq(dt);
    double tmp = 1. + sq(Wc_dt[n]) / 4.;
    w02 += wp_dt2[n] / tmp;
    gam2 += wp_dt2[n] * sq(Wc_dt[n]) / tmp;
    delta += wp_dt2[n] * Wc_dt[n] / tmp;
    Eigen::Vector3d J(f[MX] * qbym[n], f[MY] * qbym[n], f[MZ] * qbym[n]);
    K -= dt / tmp * (J + sq(Wc_dt[n] / 2.) * b * b.dot(J) - (Wc_dt[n] / 2.) * b.cross(J));
  }
  double Delta2 = sq(delta) / (1. + w02 / 4.);

  Eigen::Vector3d F(em[EX] * epsilon0, em[EY] * epsilon0, em[EZ] * epsilon0);
  Eigen::Vector3d F_halfK = F + 0.5 * K;

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

  double chargeDens = 0.0;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    double keOld = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
    Eigen::Vector3d J(f[MX] * qbym[n], f[MY] * qbym[n], f[MZ] * qbym[n]);
    Eigen::Vector3d Jstar = J + Fbar * (wp_dt2[n] / dt / 2.);
    Eigen::Vector3d J_new = 2. * (Jstar + sq(Wc_dt[n] / 2.) * b * b.dot(Jstar) - (Wc_dt[n] / 2.) * b.cross(Jstar)) / (1. + sq(Wc_dt[n] / 2.)) - J;

    f[MX] = J_new[0] / qbym[n];
    f[MY] = J_new[1] / qbym[n];
    f[MZ] = J_new[2] / qbym[n];
    if (sd->hasPressure)
    {
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld;
    }
    chargeDens += qbym[n] * f[RHO];
  } 

  //------------> update correction potential
  double crhoc = sd->chi_e * chargeDens/sd->epsilon0;
  em[PHIE] += dt * crhoc;
}
