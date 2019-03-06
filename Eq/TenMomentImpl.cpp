// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment core functions
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cmath>
#include <TenMomentImpl.h>

static double sq(double x) { return x*x; }

/* Convert conserved to primitive variables */
static void primitive(const double q[10], double out[10]) 
{ 
   out[0] = q[0]; 
   out[1] = q[1]/q[0]; 
   out[2] = q[2]/q[0]; 
   out[3] = q[3]/q[0]; 
   out[4] = q[4]-sq(q[1])/q[0]; 
   out[5] = q[5]-(q[1]*q[2])/q[0]; 
   out[6] = q[6]-(q[1]*q[3])/q[0]; 
   out[7] = q[7]-sq(q[2])/q[0]; 
   out[8] = q[8]-(q[2]*q[3])/q[0]; 
   out[9] = q[9]-sq(q[3])/q[0]; 
} 

/* Multiply by phi prime */
static void mulByPhiPrime(double p0, double u1, double u2, double u3, const double w[10], double out[10]) 
{ 
   out[0] = w[0]; 
   out[1] = w[0]*u1+w[1]*p0; 
   out[2] = w[0]*u2+w[2]*p0; 
   out[3] = w[0]*u3+w[3]*p0; 
   out[4] = w[0]*sq(u1)+2*w[1]*p0*u1+w[4]; 
   out[5] = w[0]*u1*u2+w[1]*p0*u2+w[2]*p0*u1+w[5]; 
   out[6] = w[0]*u1*u3+w[1]*p0*u3+w[3]*p0*u1+w[6]; 
   out[7] = w[0]*sq(u2)+2*w[2]*p0*u2+w[7]; 
   out[8] = w[0]*u2*u3+w[2]*p0*u3+w[3]*p0*u2+w[8]; 
   out[9] = w[0]*sq(u3)+2*w[3]*p0*u3+w[9]; 
} 

// Riemann problem for Ten-moment equations: `delta` is the vector we
// wish to split, `ql`/`qr` the left/right states. On output, `waves`
// and `s` contain the waves and speeds. waves is a mwave X meqn
// matrix. See LeVeque's book for explanations.  Note: This code is
// essentially based on code used in my (Ammar) thesis i.e.
// Miniwarpx. Parts of the code were generated automatically from
// Maxima.
void gkylTenMomentRp(int dir, double delta[10], double ql[10], double qr[10], double *waves, double s[5])
{
  // build permutation arrays depending on RP direction. Pressure
  // tensor permutation is built by permuting (i,j) indices in the
  // same way as the i and j velocities. These are 1-indexed as I am
  // already horribly confused with 0- and 1-based indexing.
  int d[4], dp[7];
  if (dir == 1)
  {
    d[1] = 1; d[2] = 2; d[3] = 3;
    dp[1] = 4; dp[2] = 5; dp[3] = 6; dp[4] = 7; dp[5] = 8; dp[6] = 9;
  }
  else if (dir == 2)
  {
    d[1] = 2; d[2] = 3; d[3] = 1;
    dp[1] = 7; dp[2] = 8; dp[3] = 5; dp[4] = 9; dp[5] = 6; dp[6] = 4;
  }
  else if (dir == 3)
  {
    d[1] = 3; d[2] = 1; d[3] = 2;
    dp[1] = 9; dp[2] = 6; dp[3] = 8; dp[4] = 4; dp[5] = 5; dp[6] = 7;
  }

  double vl[10], vr[10];
  primitive(ql, vl);
  primitive(qr, vr);

  // compute Roe averages
  double sqrl = std::sqrt(vl[0]), sqrr = std::sqrt(vr[0]);
  double sqr1 = 1/(sqrl+sqrr);
  
  double p0 = sqrl*sqrr;
  double p2s1 = sq(p0*sqr1);
  
  double u1 = (sqrl*vl[d[1]] + sqrr*vr[d[1]])*sqr1;
  double u2 = (sqrl*vl[d[2]] + sqrr*vr[d[2]])*sqr1;
  double u3 = (sqrl*vl[d[3]] + sqrr*vr[d[3]])*sqr1;
  double p11 = (sqrr*vl[dp[1]]+sqrl*vr[dp[1]])*sqr1 + 1.0/3.0*p2s1*(vr[d[1]]-vl[d[1]])*(vr[d[1]]-vl[d[1]]);
  double p12 = (sqrr*vl[dp[2]]+sqrl*vr[dp[2]])*sqr1 + 1.0/3.0*p2s1*(vr[d[1]]-vl[d[1]])*(vr[d[2]]-vl[d[2]]);
  double p13 = (sqrr*vl[dp[3]]+sqrl*vr[dp[3]])*sqr1 + 1.0/3.0*p2s1*(vr[d[1]]-vl[d[1]])*(vr[d[3]]-vl[d[3]]);
  double p22 = (sqrr*vl[dp[4]]+sqrl*vr[dp[4]])*sqr1 + 1.0/3.0*p2s1*(vr[d[2]]-vl[d[2]])*(vr[d[2]]-vl[d[2]]);
  double p23 = (sqrr*vl[dp[5]]+sqrl*vr[dp[5]])*sqr1 + 1.0/3.0*p2s1*(vr[d[2]]-vl[d[2]])*(vr[d[3]]-vl[d[3]]);
  double p33 = (sqrr*vl[dp[6]]+sqrl*vr[dp[6]])*sqr1 + 1.0/3.0*p2s1*(vr[d[3]]-vl[d[3]])*(vr[d[3]]-vl[d[3]]);

  // for multiplication by phi' we need to use unrotated values
  double v[4];
  v[d[1]] = u1; v[d[2]] = u2; v[d[3]] = u3;

  double phiDelta[10];

  // pre-multiply jump (delta) by phiPrime inverse: we do this as
  // jumps are in conserved variables, while left eigenvectors used
  // below are computed from primitive variables
  phiDelta[0] = delta[0];
  phiDelta[1] = delta[d[1]]/p0-(1.0*delta[0]*u1)/p0; 
  phiDelta[2] = delta[d[2]]/p0-(1.0*delta[0]*u2)/p0; 
  phiDelta[3] = delta[d[3]]/p0-(1.0*delta[0]*u3)/p0; 
  phiDelta[4] = delta[0]*sq(u1)-2.0*delta[d[1]]*u1+delta[dp[1]]; 
  phiDelta[5] = delta[0]*u1*u2-1.0*delta[d[1]]*u2-1.0*delta[d[2]]*u1+delta[dp[2]]; 
  phiDelta[6] = delta[0]*u1*u3-1.0*delta[d[1]]*u3-1.0*delta[d[3]]*u1+delta[dp[3]]; 
  phiDelta[7] = delta[0]*sq(u2)-2.0*delta[d[2]]*u2+delta[dp[4]]; 
  phiDelta[8] = delta[0]*u2*u3-1.0*delta[d[2]]*u3-1.0*delta[d[3]]*u2+delta[dp[5]]; 
  phiDelta[9] = delta[0]*sq(u3)-2.0*delta[d[3]]*u3+delta[dp[6]];

  // predefine some constants
  double p11sq = sq(p11), p12sq = sq(p12), p13sq = sq(p13), p11th = std::pow(p11, 1.5);
  double sqp0 = std::sqrt(p0), sqp11 = std::sqrt(p11);
  
  double leftProj[10];
  // project jumps on left eigenvectors (pray that C++ compiler eliminates common subexpressions) [Gen from Maxima]
  leftProj[0] = (0.5*phiDelta[1]*sqp0*sqp11*p12)/p11sq-(0.5*phiDelta[4]*p12)/p11sq-(0.5*phiDelta[2]*sqp0)/sqp11+(0.5*phiDelta[5])/p11; 
  leftProj[1] = (0.5*phiDelta[1]*sqp0*p13)/p11th-(0.5*phiDelta[4]*p13)/p11sq-(0.5*phiDelta[3]*sqp0)/sqp11+(0.5*phiDelta[6])/p11; 
  leftProj[2] = (-(0.5*phiDelta[1]*sqp0*sqp11*p12)/p11sq)-(0.5*phiDelta[4]*p12)/p11sq+(0.5*phiDelta[2]*sqp0)/sqp11+(0.5*phiDelta[5])/p11; 
  leftProj[3] = (-(0.5*phiDelta[1]*sqp0*p13)/p11th)-(0.5*phiDelta[4]*p13)/p11sq+(0.5*phiDelta[3]*sqp0)/sqp11+(0.5*phiDelta[6])/p11; 
  leftProj[4] = (0.16666666666666666667*phiDelta[4])/p11sq-(0.2886751345948129*phiDelta[1]*sqp0)/p11th; 
  leftProj[5] = (0.2886751345948129*phiDelta[1]*sqp0)/p11th+(0.16666666666666666667*phiDelta[4])/p11sq; 
  leftProj[6] = phiDelta[0]-(0.3333333333333333333*phiDelta[4]*p0)/p11; 
  leftProj[7] = (-(0.3333333333333333333*phiDelta[4]*p11*p22)/p11sq)+(1.333333333333333333*phiDelta[4]*p12sq)/p11sq-(2.0*phiDelta[5]*p12)/p11+phiDelta[7]; 
  leftProj[8] = (-(0.3333333333333333333*phiDelta[4]*p11*p23)/p11sq)+(1.333333333333333333*phiDelta[4]*p12*p13)/p11sq-(1.0*phiDelta[5]*p13)/p11-(1.0*phiDelta[6]*p12)/p11+phiDelta[8]; 
  leftProj[9] = (-(0.3333333333333333333*phiDelta[4]*p11*p33)/p11sq)+(1.333333333333333333*phiDelta[4]*p13sq)/p11sq-(2.0*phiDelta[6]*p13)/p11+phiDelta[9]; 

  // compute waves and speeds
  double wv[10];

  // Wave 1: (ev 1 and 2 are repeated)
  s[0] = u1-std::sqrt(p11/p0);
  wv[0] = 0.0; 
  wv[d[1]] = 0.0; 
  wv[d[2]] = -(1.0*leftProj[0]*sqp11)/sqp0; 
  wv[d[3]] = -(1.0*leftProj[1]*sqp11)/sqp0; 
  wv[dp[1]] = 0.0; 
  wv[dp[2]] = leftProj[0]*p11; 
  wv[dp[3]] = leftProj[1]*p11; 
  wv[dp[4]] = 2.0*leftProj[0]*p12; 
  wv[dp[5]] = leftProj[0]*p13+leftProj[1]*p12; 
  wv[dp[6]] = 2.0*leftProj[1]*p13;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[0]);

  // Wave 2: (ev 3 and 4 are repeated)
  s[1] = u1+std::sqrt(p11/p0);
  wv[0] = 0.0; 
  wv[d[1]] = 0.0; 
  wv[d[2]] = (leftProj[2]*sqp11)/sqp0; 
  wv[d[3]] = (leftProj[3]*sqp11)/sqp0; 
  wv[dp[1]] = 0.0; 
  wv[dp[2]] = leftProj[2]*p11; 
  wv[dp[3]] = leftProj[3]*p11; 
  wv[dp[4]] = 2.0*leftProj[2]*p12; 
  wv[dp[5]] = leftProj[2]*p13+leftProj[3]*p12; 
  wv[dp[6]] = 2.0*leftProj[3]*p13;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[10]);
  
  // Wave 3 (ev 5)
  s[2] = u1-std::sqrt(3*p11/p0);
  wv[0] = leftProj[4]*p0*p11; 
  wv[d[1]] = -(1.732050807568877*leftProj[4]*p11th)/sqp0; 
  wv[d[2]] = -(1.732050807568877*leftProj[4]*sqp11*p12)/sqp0; 
  wv[d[3]] = -(1.732050807568877*leftProj[4]*sqp11*p13)/sqp0; 
  wv[dp[1]] = 3.0*leftProj[4]*p11sq; 
  wv[dp[2]] = 3.0*leftProj[4]*p11*p12; 
  wv[dp[3]] = 3.0*leftProj[4]*p11*p13; 
  wv[dp[4]] = leftProj[4]*p11*p22+2.0*leftProj[4]*p12sq; 
  wv[dp[5]] = leftProj[4]*p11*p23+2.0*leftProj[4]*p12*p13; 
  wv[dp[6]] = leftProj[4]*p11*p33+2.0*leftProj[4]*p13sq;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[20]);

  // Wave 4 (ev 6)
  s[3] = u1+std::sqrt(3*p11/p0);
  wv[0] = leftProj[5]*p0*p11; 
  wv[d[1]] = (1.732050807568877*leftProj[5]*p11th)/sqp0; 
  wv[d[2]] = (1.732050807568877*leftProj[5]*sqp11*p12)/sqp0; 
  wv[d[3]] = (1.732050807568877*leftProj[5]*sqp11*p13)/sqp0; 
  wv[dp[1]] = 3.0*leftProj[5]*p11sq; 
  wv[dp[2]] = 3.0*leftProj[5]*p11*p12; 
  wv[dp[3]] = 3.0*leftProj[5]*p11*p13; 
  wv[dp[4]] = leftProj[5]*p11*p22+2.0*leftProj[5]*p12sq; 
  wv[dp[5]] = leftProj[5]*p11*p23+2.0*leftProj[5]*p12*p13; 
  wv[dp[6]] = leftProj[5]*p11*p33+2.0*leftProj[5]*p13sq;

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[30]);

  // Wave 5: (ev 7, 8, 9, 10 are repeated)
  s[4] = u1;
  wv[0] = leftProj[6]; 
  wv[d[1]] = 0.0; 
  wv[d[2]] = 0.0; 
  wv[d[3]] = 0.0; 
  wv[dp[1]] = 0.0; 
  wv[dp[2]] = 0.0; 
  wv[dp[3]] = 0.0; 
  wv[dp[4]] = leftProj[7]; 
  wv[dp[5]] = leftProj[8]; 
  wv[dp[6]] = leftProj[9];

  mulByPhiPrime(p0, v[1], v[2], v[3], wv, &waves[40]);
}

// Qfluctuations for ten-moment equations
void gkylTenMomentQFluct(double *waves, double s[5], double amdq[10], double apdq[10])
{
  for (unsigned m=0; m<10; ++m)
  {
    amdq[m] = apdq[m] = 0.0;
    for (unsigned mw=0; mw<5; ++mw)
    {
      if (s[mw] < 0)
        amdq[m] += s[mw]*waves[m+10*mw];
      else
        apdq[m] += s[mw]*waves[m+10*mw];
    }
  }
}
