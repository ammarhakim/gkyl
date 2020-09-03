// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for use in multimomlinear Tool. This code computes
// eigenvalues and eigenvectors of the complex matrix supplied to it.
// 
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// std includes
#include <complex>

// eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace Gkyl { class EigenEigen; }

extern "C" {
    // C wrappers for interfacing with EigenEigen class
    void* new_EigenEigen(int N, double *mre, double *mim);
    void delete_EigenEigen(Gkyl::EigenEigen* ee);
    void compute_EigenEigen(Gkyl::EigenEigen* ee, double *ere, double *eim);
}

namespace Gkyl {
  
  // Computes eigenvalues of supplied complex matrix
  class EigenEigen {
    public:
      /**
       * Construct a new object to compute eigenvalues of given
       * complex matrix.
       *
       * @param N Rows and columns of square matrix
       * @param mre Real part of matrix, arranged in row-major order
       * @param mim Imaginary part of matrix, arranged in row-major order
       */
      EigenEigen(int N, double *mre, double *mim);

      /**
       * Computes eigenvalues of system.
       *
       * @param ere On output, stores real part of eigenvalues
       * @param eim On outout, stores imaginary part of eigenvalues
       */
      void compute(double *ere, double *eim);
      
    private:
      /** Dispersion matrix (complex double type) */
      Eigen::MatrixXcd D;
  };
}
