extern "C" {
  void reflectionLoop(double* fItr, double* ordinates, double* weights, double* iWeights, double* iOrdinates, double* negBasis, double* posBasis, double* vc, double* dv, double C, double elemCharge, double me, double electronAffinity, double effectiveMass, int vdim, int cdim, int numOrd, int N, int numBasis, int ignore);
  double getR(double E, double xi, double electronAffinity, double effectiveMass, double C, int N, double* iWeights, double* iOrdinates);
  double quad(double* f, double a, double b, double* weights, int numOrd);
  double getT(double E, double xi, double electronAffinity, double effectiveMass);
  double getE(double* v, double* vc, double* dv, double me, double elemCharge, int vdim);
  double getXi(double* v, double* vc, double* dv, int vdim);
  double getXic(double E, double electronAffinity, double effectiveMass);
};
