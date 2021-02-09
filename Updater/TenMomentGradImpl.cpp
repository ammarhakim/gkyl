// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment gradient-based closure
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <TenMomentGradImpl.h>
#include <cmath>
#include <stdio.h>

void
gkylTenMomentHeatFlux(const double alpha, const double* dT1, const double* dT2, const double* dT3, const double* f, double* q)
{
  double r = f[0];
  double pxx = f[4];
  double pyy = f[7];
  double pzz = f[9];

  double p = (pxx+pyy+pzz)/3.0;
  double vt = std::sqrt(p/r);

  q[Q111] = alpha*vt*r*(dT1[T11] + dT1[T11] + dT1[T11])/3.0;
  q[Q112] = alpha*vt*r*(dT1[T12] + dT1[T12] + dT2[T11])/3.0;
  q[Q113] = alpha*vt*r*(dT1[T13] + dT1[T13] + dT3[T11])/3.0;
  q[Q122] = alpha*vt*r*(dT1[T22] + dT2[T12] + dT2[T12])/3.0;
  q[Q123] = alpha*vt*r*(dT1[T23] + dT2[T13] + dT3[T12])/3.0;
  q[Q133] = alpha*vt*r*(dT1[T33] + dT3[T13] + dT3[T13])/3.0; 
  q[Q222] = alpha*vt*r*(dT2[T22] + dT2[T22] + dT2[T22])/3.0;
  q[Q223] = alpha*vt*r*(dT2[T23] + dT2[T23] + dT3[T22])/3.0;
  q[Q233] = alpha*vt*r*(dT2[T33] + dT3[T23] + dT3[T23])/3.0;
  q[Q333] = alpha*vt*r*(dT3[T33] + dT3[T33] + dT3[T33])/3.0;
}

void 
gkylTenMomentAccumulateGradClosure(const double* divQ1, const double* divQ2, const double* divQ3, double* f)
{
  // compute updated pressure tensor component
  f[4] = (divQ1[0] + divQ2[0] + divQ3[0]);
  f[5] = (divQ1[1] + divQ2[1] + divQ3[1]);
  f[6] = (divQ1[2] + divQ2[2] + divQ3[2]);
  f[7] = (divQ1[3] + divQ2[3] + divQ3[3]);
  f[8] = (divQ1[4] + divQ2[4] + divQ3[4]);
  f[9] = (divQ1[5] + divQ2[5] + divQ3[5]);
}

void
gkylTenMomentGradT(const int dir, const double* dxv, const double* fL, const double* fR, double* dT)
{
  // fetch grid spacing in given direction
  const double dx = dxv[dir];
  double rL = fL[0];
  double TxxL = fL[4]/rL;
  double TxyL = fL[5]/rL;
  double TxzL = fL[6]/rL;
  double TyyL = fL[7]/rL;
  double TyzL = fL[8]/rL;
  double TzzL = fL[9]/rL;

  double rR = fR[0];
  double TxxR = fR[4]/rR;
  double TxyR = fR[5]/rR;
  double TxzR = fR[6]/rR;
  double TyyR = fR[7]/rR;
  double TyzR = fR[8]/rR;
  double TzzR = fR[9]/rR;

  dT[0] = 0.5*(TxxR - TxxL)/dx;
  dT[1] = 0.5*(TxyR - TxyL)/dx;
  dT[2] = 0.5*(TxzR - TxzL)/dx;
  dT[3] = 0.5*(TyyR - TyyL)/dx;
  dT[4] = 0.5*(TyzR - TyzL)/dx;
  dT[5] = 0.5*(TzzR - TzzL)/dx;
}

void
gkylTenMomentDivQX(const double* dxv, const double* qL, const double* qR, double* divQ)
{
  // fetch grid spacing in x direction
  const double dx = dxv[0];
  divQ[0] = 0.5*(qR[Q111] - qL[Q111])/dx;
  divQ[1] = 0.5*(qR[Q112] - qL[Q112])/dx;
  divQ[2] = 0.5*(qR[Q113] - qL[Q113])/dx;
  divQ[3] = 0.5*(qR[Q122] - qL[Q122])/dx;
  divQ[4] = 0.5*(qR[Q123] - qL[Q123])/dx;
  divQ[5] = 0.5*(qR[Q133] - qL[Q133])/dx;
}

void
gkylTenMomentDivQY(const double* dxv, const double* qL, const double* qR, double* divQ)
{
  // fetch grid spacing in y direction
  const double dx = dxv[1];
  divQ[0] = 0.5*(qR[Q112] - qL[Q112])/dx;
  divQ[1] = 0.5*(qR[Q122] - qL[Q122])/dx;
  divQ[2] = 0.5*(qR[Q123] - qL[Q123])/dx;
  divQ[3] = 0.5*(qR[Q222] - qL[Q222])/dx;
  divQ[4] = 0.5*(qR[Q223] - qL[Q223])/dx;
  divQ[5] = 0.5*(qR[Q233] - qL[Q233])/dx;
}

void
gkylTenMomentDivQZ(const double* dxv, const double* qL, const double* qR, double* divQ)
{
  // fetch grid spacing in z direction
  const double dx = dxv[2];
  divQ[0] = 0.5*(qR[Q113] - qL[Q113])/dx;
  divQ[1] = 0.5*(qR[Q123] - qL[Q123])/dx;
  divQ[2] = 0.5*(qR[Q133] - qL[Q133])/dx;
  divQ[3] = 0.5*(qR[Q223] - qL[Q223])/dx;
  divQ[4] = 0.5*(qR[Q233] - qL[Q233])/dx;
  divQ[5] = 0.5*(qR[Q333] - qL[Q333])/dx;
}