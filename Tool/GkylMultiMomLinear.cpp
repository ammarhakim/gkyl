// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for use in multimomlinear Tool. This code computes
// eigenvalues and eigenvectors of the complex matrix supplied to it.
// 
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylMultiMomLinear.h>
#include <iostream>

namespace Gkyl {
  EigenEigen::EigenEigen(int N, double *mre, double *mim, int cv)
    : D(N,N), calcVec(cv == 0 ? false : true) {
    for (unsigned r=0; r<N; ++r)
      for (unsigned c=0; c<N; ++c) {
        D(r,c).real(mre[r*N+c]);
        D(r,c).imag(mim[r*N+c]);
      }
  }

  void EigenEigen::compute(double *ere, double *eim) {
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;

    // compute the eigensystem
    solver.compute(D, calcVec);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd>::EigenvalueType ev
      = solver.eigenvalues();

    Eigen::ComputationInfo info = solver.info();
    if (info != Eigen::ComputationInfo::Success) {
      // should throw an error
      std::cout << "!!!***!!! Eigenvalue computation failed!" << std::endl;
    }

    // copy into output arrays (ev is a column vector of eigenvalues)
    for (unsigned i=0; i<ev.rows(); ++i) {
      ere[i] = ev(i,0).real();
      eim[i] = ev(i,0).imag();
    }
  }
  
}

extern "C" void* new_EigenEigen(int N, double *mre, double *mim, int calcVec) {
  Gkyl::EigenEigen *ee = new Gkyl::EigenEigen(N, mre, mim, calcVec);
  return reinterpret_cast<void*>(ee);
}

extern "C" void delete_EigenEigen(Gkyl::EigenEigen* ee) {
  delete ee;
}

extern "C" void compute_EigenEigen(Gkyl::EigenEigen* ee, double *ere, double *eim) {
  ee->compute(ere, eim);
}


