// Gkyl ------------------------------------------------------------------------
//
// Euler RP in C++
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <GkEulerEqn.h>

// Makes indexing easier
#define IDX(arr, i, j) arr[5*i+j]

static inline double pressure(const EulerEqn_t *e, const double *q) {
  return (e->gasGamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
}

void
Euler_rp(EulerEqn_t *e, int dir, double *delta, double *ql, double *qr, double *waves, double *s)
{
  // set up reshuffle indices for various direction Riemann problem
  int d[3];
  if (dir == 1) {
    d[0] = 1; d[1] = 2; d[2] = 3;
  }
  else if (dir == 2) {
    d[0] = 2; d[1] = 3; d[2] = 1;
  }
  else {
    d[0] = 3; d[1] = 1; d[2] = 2;
  }

  double g1 = e->gasGamma-1;
  double rhol = ql[0], rhor = qr[0];
  double pl = pressure(e,ql), pr = pressure(e,qr);

  // Roe averages: see Roe's original 1986 paper or LeVeque book
  double srrhol = std::sqrt(rhol), srrhor = std::sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[d[0]]*ravgl1 + qr[d[0]]*ravgr1)*ravg2;
  double v = (ql[d[1]]*ravgl1 + qr[d[1]]*ravgr1)*ravg2;
  double w = (ql[d[2]]*ravgl1 + qr[d[2]]*ravgr1)*ravg2;
  double enth = ((ql[4]+pl)*ravgl1 + (qr[4]+pr)*ravgr1)*ravg2;

  // See http://ammar-hakim.org/sj/euler-eigensystem.html for notation
  // and meaning of these terms
  double q2 = u*u+v*v+w*w;
  double aa2 = g1*(enth-0.5*q2);
  double a = std::sqrt(aa2);
  double g1a2 = g1/aa2, euv = enth-q2;
  
  // compute projections of jump
  double a4 = g1a2*(euv*delta[0] + u*delta[d[0]] + v*delta[d[1]] + w*delta[d[2]] - delta[4]);
  double a2 = delta[d[1]] - v*delta[0];
  double a3 = delta[d[2]] - w*delta[0];
  double a5 = 0.5*(delta[d[0]] + (a-u)*delta[0] - a*a4)/a;
  double a1 = delta[0] - a4 - a5;

  // wave 1: eigenvalue is u-c
  double *wv = &waves[5*0];
  wv[0] = a1;
  wv[d[0]] = a1*(u-a);
  wv[d[1]] = a1*v;
  wv[d[2]] = a1*w;
  wv[4] = a1*(enth-u*a);
  s[0] = u-a;

  // wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5*1];
  wv[0] = a4;
  wv[d[0]] = a4*u;
  wv[d[1]] = a4*v + a2;
  wv[d[2]] = a4*w + a3;
  wv[4] = a4*0.5*q2 + a2*v + a3*w;
  s[1] = u;

  // wave 3: eigenvalue is u+c
  wv = &waves[5*2];
  wv[0] = a5;
  wv[d[0]] = a5*(u+a);
  wv[d[1]] = a5*v;
  wv[d[2]] = a5*w;
  wv[4] = a5*(enth+u*a);
  s[2] = u+a;
}
