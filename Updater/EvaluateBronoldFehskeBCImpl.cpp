#include <EvaluateBronoldFehskeBCImpl.h>
#include <math.h>

void reflectionLoop(double* fItr, double* ordinates, double* weights, double* iWeights, double* iOrdinates, double* negBasis, double* posBasis, double* vc, double* dv, double C, double elemCharge, double me, double electronAffinity, double effectiveMass, int vdim, int cdim, int numOrd, int N, int numBasis, int ignore)
{
  if((ignore == -1 && vc[0] < 0.0) || (ignore == 1 && vc[0] > 0.0)) {
    for(int m=0; m<numBasis; m++) {
      for(int n=0; n<numBasis; n++) {
        fItr[numBasis*m + n] = 0.0;
      }
    }
  }
  else {
    for(int m=0; m<numBasis; m++) {
      for(int n=0; n<numBasis; n++) {
        double Rquad[numOrd] = {};
        for(int i=0; i<numOrd; i++) {
          double v[vdim] = {};
	  for(int d=0; d<vdim; d++) {
	    v[d] = ordinates[i*(vdim + cdim) + (d + cdim)];
	  }
	  double E = getE(v, vc, dv, me, elemCharge, vdim);
	  double xi = getXi(v, vc, dv, vdim);
	  double nB = negBasis[i*numBasis + n];
	  double pB = posBasis[i*numBasis + m];
	  Rquad[i] = getR(E, xi, electronAffinity, effectiveMass, C, N, iWeights, iOrdinates)*nB*pB;
        }
        fItr[numBasis*m + n] = quad(Rquad, -1.0, 1.0, weights, numOrd);
      }
    }
  }
}

double getR(double E, double xi, double electronAffinity, double effectiveMass, double C, int N, double* iWeights, double* iOrdinates)
{
  double xic = getXic(E, electronAffinity, effectiveMass);
  double Tr = getT(E, xi, electronAffinity, effectiveMass);
  double Tquad [N] = {};
  for(int i=0; i<N; i++) {
    Tquad[i] = getT(E, (1 - xic)*iOrdinates[i]/2 + (1 + xic)/2, electronAffinity, effectiveMass);
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
  }
  return ans;
}

double getT(double E, double xi, double electronAffinity, double effectiveMass)
{
  double eta = sqrt(1 - ((E - electronAffinity)/(effectiveMass*E))*(1 - pow(xi, 2)));
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

double getXic(double E, double electronAffinity, double effectiveMass)
{
  if(E < electronAffinity/(1 - effectiveMass)) {
    return 0.0;
  }
  else {
    return sqrt(1 - ((effectiveMass*E)/(E - electronAffinity)));
  }
}
