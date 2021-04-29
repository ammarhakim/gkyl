#include <EvaluateBronoldFehskeBCImpl.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

void reflectionLoop(double* fItr, double* ordinates, double* weights, double* iWeights, double* negBasis, double* posBasis, double* vc, double* dv, double C, double elemCharge, double me, double electronAffinity, double effectiveMass, int vdim, int cdim, int numOrd, int N, int numBasis)
{
  for(int m=0; m<numBasis; m++) {
    for(int n=0; n<numBasis; n++) {
      double Rquad[numOrd] = {};
      for(int i=0; i<numOrd; i++) {
	double v[vdim] = {};
	for(int d=0; d<vdim; d++) {
	  v[d] = ordinates[i*numOrd + (d + cdim)];
	}
	double E = getE(v, vc, dv, me, elemCharge, vdim);
	double xi = getXi(v, vc, dv, vdim);
	double nB = negBasis[i*numOrd + n];
	double pB = posBasis[i*numOrd + m];
	Rquad[i] = getR(E, xi, electronAffinity, effectiveMass, C, N, iWeights)*nB*pB;
	//std::cout << i << "\n";
      }
      fItr[numBasis*(m - 1) + n] = quad(Rquad, -1.0, 1.0, weights, numOrd);
      std::cout << numBasis*(m - 1) + n << ": " << fItr[numBasis*(m - 1) + n] << "\n";
    }
  }
}

double getR(double E, double xi, double electronAffinity, double effectiveMass, double C, int N, double* iWeights)
{
  double xic;
  if(E < electronAffinity/(1 - effectiveMass)) {
    xic = 0.0;
  }
  else {
    xic = sqrt(1 - ((effectiveMass*E)/(E - electronAffinity)));
  }
  double Tr = getT(E, xi, electronAffinity, effectiveMass);
  double Tquad [N] = {};
  for(int i=0; i<N; i++) {
    Tquad[i] = getT(E, xi, electronAffinity, effectiveMass);
  }
  double Tint = quad(Tquad, xic, 1.0, iWeights, N);
  if(E < electronAffinity) {
    return 1.0;
  }
  else {
    if(xi > xic) {
      return 1 - (Tr/(1 + (C/xi))) - ((C/xi)/(1 + (C/xi)))*Tint;
    }
    else {
      return 1 - ((C/xi)/(1 + (C/xi)))*Tint;
    }
  }
}

double quad(double* f, double a, double b, double* weights, int numOrd)
{
  double ans = 0;
  for(int i=0; i<numOrd; i++) {
    ans += (b - a)*(weights[i]*f[i]/2);
    //std::cout << "C++ [f, weights, ans]: " << f[i] << "  " << weights[i] << "  " << ans << "\n";
  }
  return ans;
}

double getT(double E, double xi, double electronAffinity, double effectiveMass)
{
  double eta = sqrt(1 - ((E - electronAffinity)/(effectiveMass*E))*(1 - pow(xi, 2)));
  //std::cout << sqrt(((E - electronAffinity)/(effectiveMass))*(1 - pow(xi, 2))) << "\n";
  return (4*effectiveMass*sqrt(E - electronAffinity)*xi*sqrt(effectiveMass*E)*eta)/pow((effectiveMass*sqrt(E - electronAffinity)*xi + sqrt(effectiveMass*E)*eta), 2);
}

double getE(double v[], double* vc, double* dv, double me, double elemCharge, int vdim)
{
  double ans = 0;
  for(int d=0; d<vdim; d++) {
    ans += 0.5*me*pow(vc[d] + v[d]*(dv[d]/2), 2)/elemCharge;
  }
  return ans;
}

double getXi(double v[], double* vc, double* dv, int vdim)
{
  double ans = 0;
  for(int d=0; d<vdim; d++) {
    ans += pow(vc[d] + v[d]*(dv[d]/2), 2);
  }
  return abs(vc[0] + v[0]*(dv[0]/2))/sqrt(ans);
}
