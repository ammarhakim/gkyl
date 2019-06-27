#include <EvalOnNodesImpl.h>
#include <FemMatrices.h>

using namespace Eigen;

void nodToMod(double* fN, int numNodes, int numVal, int ndim, int p, double* fM)
{
  Map<MatrixXd, RowMajor> fNod(fN, numNodes, numVal);
  Map<MatrixXd, RowMajor> fMod(fM, numNodes, numVal);

  MatrixXd nodToMod = MatrixXd::Zero(numNodes, numNodes);
  getNodToModMatrix(nodToMod, ndim, p);

  fMod = nodToMod*fNod;
}

