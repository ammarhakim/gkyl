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
static double cube(double x) { return (x)*(x)*(x); }
static double sgn(double x) {
  if (x>0) return 1;
  else if (x < 0) return -1;
  else return 0;
}

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

    if (fd[n].evolve)
    {
      // off-diagonal elements of lhs
      // eqn. for X-component of current
      lhs(fidx(n,X), fidx(n,Y)) = -dt1*qbym*(em[BZ]+staticEm[BZ]);
      lhs(fidx(n,X), fidx(n,Z)) = dt1*qbym*(em[BY]+staticEm[BY]);
      lhs(fidx(n,X), eidx(X)) = -dt1*qbym2*f[RHO];

      // eqn. for Y-component of current
      lhs(fidx(n,Y), fidx(n,X)) = dt1*qbym*(em[BZ]+staticEm[BZ]);
      lhs(fidx(n,Y), fidx(n,Z)) = -dt1*qbym*(em[BX]+staticEm[BX]);
      lhs(fidx(n,Y), eidx(Y)) = -dt1*qbym2*f[RHO];

      // eqn. for Z-component of current
      lhs(fidx(n,Z), fidx(n,X)) = -dt1*qbym*(em[BY]+staticEm[BY]);
      lhs(fidx(n,Z), fidx(n,Y)) = dt1*qbym*(em[BX]+staticEm[BX]);
      lhs(fidx(n,Z), eidx(Z)) = -dt1*qbym2*f[RHO];
    }

    // diagonal elements of lhs
    lhs(fidx(n,X), fidx(n,X)) = 1.0;
    lhs(fidx(n,Y), fidx(n,Y)) = 1.0;
    lhs(fidx(n,Z), fidx(n,Z)) = 1.0;

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
gkylFiveMomentSrcTimeCenteredDirect2(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
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
gkylFiveMomentSrcTimeCenteredDirect(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm)
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
  Eigen::Vector3d b(0., 0., 0.);
  if (Bmag > 0.)
  {
    b[0] = Bx / Bmag;
    b[1] = By / Bmag;
    b[2] = Bz / Bmag;
  }

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

  double chargeDens = 0.0;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    chargeDens += qbym[n] * f[RHO];
    if (!fd[n].evolve)
      continue;
    double keOld = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
    Eigen::Vector3d Jstar = J[n] + Fbar * (wp_dt2[n] / dt / 2.);
    Eigen::Vector3d J_new = 2. * (Jstar + sq(Wc_dt[n] / 2.) * b * b.dot(Jstar) - (Wc_dt[n] / 2.) * b.cross(Jstar)) / (1. + sq(Wc_dt[n] / 2.)) - J[n];

    f[MX] = J_new[0] / qbym[n];
    f[MY] = J_new[1] / qbym[n];
    f[MZ] = J_new[2] / qbym[n];
    if (sd->hasPressure)
    {
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld;
    }
  } 

  //------------> update correction potential
  double crhoc = sd->chi_e * chargeDens/sd->epsilon0;
  em[PHIE] += dt * crhoc;
}

/* Updater for the parallel component of the exact source term.
 *
 * @param q_par Normalized parallel electric field and current along the
 *    background B field direction, i.e., [Ez, Jz0, Jz1, ..., Jzs, ...], where s
 *    denotes a species. The content of q_par will be modified in-place.
 * @param dt Time step.
 * @param wp Plasma frequency for each species.
 */
static void
update_par(Eigen::VectorXd &q_par, double dt, std::vector<double> &wp)
{
  unsigned nFluids = wp.size();
  double wp_tot = 0.;
  for (unsigned n=0; n < nFluids; ++n)
  {
    wp_tot += sq(wp[n]);
  }
  wp_tot = std::sqrt(wp_tot);

  // TODO reuse v0, v1
  // eigenvector with w=-w_p
  Eigen::VectorXd v0 = Eigen::VectorXd::Constant(nFluids + 1, 0.);
  // eigenvector with w=w_p
  Eigen::VectorXd v1 = Eigen::VectorXd::Constant(nFluids + 1, 0.);

  // eigenvectors at t=0
  v0[0] = wp_tot;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    v1[n + 1] = wp[n];
  }

  double coeff[2];
  coeff[0] = q_par.dot(v0) / v0.dot(v0);
  coeff[1] = q_par.dot(v1) / v1.dot(v1);

  double cost = std::cos(wp_tot * dt);
  double sint = std::sin(wp_tot * dt);

  // incremental changes to the eigenvectors
  v0[0] = wp_tot * cost - v0[0];
  v1[0] = -wp_tot * sint - v1[0];
  for (unsigned n = 0; n < nFluids; ++n)
  {
    v0[n + 1] = wp[n] * sint - v0[n + 1];
    v1[n + 1] = wp[n] * cost - v1[n + 1];
  }

  // accumulate incremental changes
  q_par += coeff[0] * v0  + coeff[1] * v1;
}


/*
 * Finding roots of a polynomial as eigenvalues of its companion matrix.
 *
 * @param coeffs The coefficients of the polynomial from the higher to lower order.
 * @return roots The roots.
 */
static Eigen::VectorXd roots(const std::vector<double> &coeffs)
{
  int N = coeffs.size() - 1;

  // companion matrix of the polynomial
  Eigen::MatrixXd M = Eigen::MatrixXd::Constant(N, N, 0);
  
  for(int n = 1; n < N; ++n){
    M(n, n-1) = 1.;
  }

  for(int n = 0; n < N; ++n){
    M(n, N-1) = -coeffs[N-n] / coeffs.front();
  }

  Eigen::VectorXd roots = M.eigenvalues().real();

  return roots;
}

/* Updater for the perpendicular component of the exact source term.
 *
 * @param q_par Normalized perpendicular electric field and current along the
 *    background B field direction, i.e., [Ex, Ey; Jx0, Jy0, Jx1, Jy1,  ...,
 *    Jxs, Jys, ...], where s denotes a species.
 * @param dt Time step.
 * @param wp Plasma frequency for each species.
 * @param Wc Signed cyclotron frequency for each species.
 *
 * @returns q_perp_. Updated vector.
 */
static Eigen::VectorXd
update_perp(Eigen::VectorXd &q_perp, double dt, std::vector<double> &wp,
            std::vector<double> &Wc)
{
  unsigned nFluids = wp.size();

  // TODO reuse v0, v1
  // compute all eigenvalues
  Eigen::VectorXd eigs(nFluids + 1);
  if (nFluids == 2)
  {
    if (false) {
      std::vector<double> poly_coeff = {
        1.,

        Wc[0] + Wc[1],

        Wc[0] * Wc[1] - sq(wp[0]) - sq(wp[1]),

        -Wc[0] * sq(wp[1]) - Wc[1] * sq(wp[0])
      };
      eigs = roots(poly_coeff);
    } else {
      // analytic solution based on
      // http://web.cs.iastate.edu/~cs577/handouts/polyroots.pdf
      double p = Wc[0] + Wc[1];
      double q = Wc[0] * Wc[1] - sq(wp[0]) - sq(wp[1]);
      double r = -Wc[0] * sq(wp[1]) - Wc[1] * sq(wp[0]);

      double a = (3 * q - sq(p)) / 3;
      double b = (2 * cube(p) - 9 * p * q + 27 * r) / 27;

      double det = sq(b) / 4 + cube(a) / 27;

      if (det < 0) {
        double tmp = 2 * std::sqrt(-a / 3);
        double phi = std::acos(-sgn(b) * std::sqrt(b * b / 4 / (-cube(a) / 27)));
        eigs[0] = tmp * std::cos((phi) / 3) - p / 3;
        eigs[1] = tmp * std::cos((phi + 2 * M_PI) / 3) - p / 3;
        eigs[2] = tmp * std::cos((phi + 4 * M_PI) / 3) - p / 3;
      } else if (det == 0) {
        double tmp = sgn(b) * std::sqrt(-a / 3);
        eigs[0] = -2 * tmp;
        eigs[1] = tmp;
        eigs[2] = tmp;
      } else {
        assert(false);
      }
    }
  } else if (nFluids == 3) {
    std::vector<double> poly_coeff = {
      1.,
      
      Wc[0] + Wc[1] + Wc[2],
      
      Wc[0] * Wc[1] + Wc[0] * Wc[2] +
      Wc[1] * Wc[2] - sq(wp[0]) - sq(wp[1]) - sq(wp[2]),

      Wc[0] * Wc[1] * Wc[2] - Wc[0] * sq(wp[1]) - Wc[0] * sq(wp[2]) -
      Wc[1] * sq(wp[0]) - Wc[1] * sq(wp[2]) - Wc[2] * sq(wp[0]) -
      Wc[2] * sq(wp[1]),
      
      -Wc[0] * Wc[1] * sq(wp[2]) -
      Wc[0] * Wc[2] * sq(wp[1]) - Wc[1] * Wc[2] * sq(wp[0])
    };
    eigs = roots(poly_coeff);
  } else {
    assert(false);
  }

  // compute the two eigenvector for each eigenvalue and accumulate their
  // contributions to the final parallel state vector
  Eigen::VectorXd q_perp_ = Eigen::VectorXd::Constant(2 * (nFluids + 1), 0.);
  Eigen::VectorXd v0(2 * (nFluids + 1));
  Eigen::VectorXd v1(2 * (nFluids + 1));
  for (unsigned n = 0; n < nFluids + 1; ++n)
  {
    double w = eigs[n];

    // compute the two eigenvectors for w at t=0
    v0[0] = 0.;
    v0[1] = 1.;
    v1[0] = 1.;
    v1[1] = 0.;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      unsigned nn = n + 1;
      double tmp = wp[n] / (w + Wc[n]);
      v0[2 * nn] = tmp;
      v0[2 * nn + 1] = 0.;
      v1[2 * nn] = 0;
      v1[2 * nn + 1] = -tmp;
    }

    // compute eigencoefficients
    double coeff[2];
    coeff[0] = q_perp.dot(v0) / v0.dot(v0);
    coeff[1] = q_perp.dot(v1) / v1.dot(v1);

    // compute the two eigenvectors for w at t=dt
    double cost = std::cos(w * dt);
    double sint = std::sin(w * dt);
    v0[0] = -sint;
    v0[1] = cost;
    v1[0] = cost;
    v1[1] = sint;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      unsigned nn = n + 1;
      double tmp = wp[n] / (w + Wc[n]);
      v0[2 * nn] = tmp * cost;
      v0[2 * nn + 1] = tmp * sint;
      v1[2 * nn] = tmp * sint;
      v1[2 * nn + 1] = -tmp * cost;
    }

    // accumulate contribution from the two eigenvectors
    q_perp_ += coeff[0] * v0  + coeff[1] * v1;
  }

  return q_perp_;
}

void
gkylFiveMomentSrcExact(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt,
                       double **ff, double *em, double *staticEm)
{
  unsigned nFluids = sd->nFluids;
  double epsilon0 = sd->epsilon0;

  std::vector<double> zeros(6, 0);
  if (!(sd->hasStatic))
  {
    staticEm = zeros.data();
  }

  Eigen::Vector3d B(
      em[BX] + staticEm[BX],
      em[BY] + staticEm[BY],
      em[BZ] + staticEm[BZ]);
  double Bmag = B.norm();

  Eigen::Vector3d E(em[EX], em[EY], em[EZ]);
  double Enorm = 1.;  // nominal normalization
  Eigen::Vector3d E_ = E / Enorm;

  std::vector<double> qbym(nFluids); 
  std::vector<double> Wc(nFluids); 
  std::vector<double> wp(nFluids);
  std::vector<double> Pnorm(nFluids); 
  std::vector<Eigen::Vector3d> J_(nFluids);
  for (unsigned n=0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    qbym[n] = fd[n].charge / fd[n].mass;
    Wc[n] = qbym[n] * Bmag;
    double wp2 = f[RHO] * sq(qbym[n]) / epsilon0;
    wp[n] = std::sqrt(wp2);

    Pnorm[n] = std::sqrt(epsilon0 * f[RHO]);
    if (fd[n].charge < 0.)
    {
      Pnorm[n] *= -1.;
    }
    J_[n][0] = f[MX] / Pnorm[n];
    J_[n][1] = f[MY] / Pnorm[n];
    J_[n][2] = f[MZ] / Pnorm[n];
  }

  // store initial kinetic energy, needed for updating final total energy
  double chargeDens = 0.0;
  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    chargeDens += qbym[n] * f[RHO];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  if (Bmag > 0.)  // TODO set threshold
  {
    double angle = std::acos(B[2] / Bmag);
    Eigen::Vector3d axis = B.cross(Eigen::Vector3d::UnitZ());
    axis.normalize();  // AngleAxisd needs a normalized axis vector
    Eigen::AngleAxisd rotate = Eigen::AngleAxisd(angle, axis);

    E_ = rotate * E_;
    for (unsigned n=0; n < nFluids; ++n)
      J_[n] = rotate * J_[n];

    // parallel component of initial condition
    Eigen::VectorXd q_par(nFluids + 1);
    q_par[0] = E_[2];
    for (unsigned n=0; n < nFluids; ++n)
    {
      q_par[n + 1] = J_[n][2];
    }

    // perpendicular component of initial condition
    Eigen::VectorXd q_perp(2 * (nFluids + 1));
    q_perp[0] = E_[0];
    q_perp[1] = E_[1];
    for (unsigned n=0; n < nFluids; ++n)
    {
      q_perp[2*(n + 1)] = J_[n][0];
      q_perp[2*(n + 1) + 1] = J_[n][1];
    }

    // update parallel component
    update_par(q_par, dt, wp);
    // update perpendicular component
    Eigen::VectorXd q_perp_ = update_perp(q_perp, dt, wp, Wc);

    // fill the full vectors
    E_[2] = q_par[0];
    for (unsigned n = 0; n < nFluids; ++n)
    {
      J_[n][2] = q_par[n + 1];
    }
    E_[0] = q_perp_[0];
    E_[1] = q_perp_[1];
    for (unsigned n = 0; n < nFluids; ++n)
    {
      unsigned nn = n + 1;
      J_[n][0] = q_perp_[2 * nn];
      J_[n][1] = q_perp_[2 * nn + 1];
    }

    // rotate back
    E_ = rotate.inverse() * E_;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      J_[n] = rotate.inverse() * J_[n];
    }

    // fill state vector
    em[EX] = E_[0] * Enorm;
    em[EY] = E_[1] * Enorm;
    em[EZ] = E_[2] * Enorm;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      double *f = ff[n];
      f[MX] = J_[n][0] * Pnorm[n];
      f[MY] = J_[n][1] * Pnorm[n];
      f[MZ] = J_[n][2] * Pnorm[n];
    }
  } else {
    Eigen::VectorXd q_par(nFluids + 1);

    for (unsigned d = 0; d < 3; ++d)
    {
      q_par[0] = E_[d];
      for (unsigned n=0; n < nFluids; ++n)
      {
        q_par[n + 1] = J_[n][d];
      }

      // FIXME avoid repeated calculation of eigenvectors
      update_par(q_par, dt, wp);

      // re-normalize back
      em[EX + d] = q_par[0] * Enorm;
      for (unsigned n = 0; n < nFluids; ++n)
      {
        double *f = ff[n];
        f[MX + d] = q_par[n + 1] * Pnorm[n];
      }
    }
  }

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 

  //------------> update correction potential
  double crhoc = sd->chi_e * chargeDens/sd->epsilon0;
  em[PHIE] += dt * crhoc;
}
