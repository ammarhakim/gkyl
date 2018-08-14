#ifndef GKYL_TEST_MapPoisson_H
#define GKYL_TEST_MapPoisson_H

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

//class creation
class MapPoisson;

extern "C"
{
  void *new_MapPoisson(int nx_, int ny_, double xl_, double yl_, double xu_, double yu_,
int xbctypel_, int ybctypel_, int xbctypeu_, int ybctypeu_, double xbcu_,
double xbcl_, double ybcu_, double ybcl_, bool sig_);
  //wrap fboth steps of solver
  void wrap_factorize(MapPoisson *d);
  void wrap_phisolve(MapPoisson *d);
  //metric wrap from lua to c
  void setMetricFuncPointer_MapPoisson(MapPoisson *d, void (*gfunc)(double *xcl, double *gl));
  //differential length functions
  void setDiffLenPointer_MapPoisson(MapPoisson *d, void (*hfunc)(double *xch, double *h));
  //mapcpp wrap from lua to c
  void setMapcpp_MapPoisson(MapPoisson *d, void (*mapcpp)(double xc, double yc, double *myxp));
  //wrap source, cond, and solution solvers
  void wrap_getSrcvalatIJ(MapPoisson *d, int i, int j, double sitrij);
  void wrap_getConvalatIJ(MapPoisson *d, int i, int j, int k, double citrij);
  double wrap_getSolvalatIJ(MapPoisson *d, int i, int j);

  //block storage struct
  typedef struct{
    Eigen::MatrixXd upper;
    Eigen::MatrixXd middle;
    Eigen::MatrixXd lower;
  } blkstor;
}


//class members
class MapPoisson {
  public:
    MapPoisson(int nx_, int ny_, double xl_, double yl_, double xu_, double yu_,
int xbctypel_, int ybctypel_, int xbctypeu_, int ybctypeu_, double xbcu_,
double xbcl_, double ybcu_, double ybcl_, bool sig_);

    //factorize laplacian
    void factorize();
    //solve system
    void phisolve();
    //metric calc pointer
    void setMetricFuncPointer(void (*gfunc)(double *xcl, double *gl)) {
      this->gfunc = gfunc;
    }
    //diff length pointer
    void setDiffLenPointer(void (*hfunc)(double *xch, double *h)) {
      this->hfunc = hfunc;
    }
    //set mapcpp function pointer
    void setMapcpp(void (*mapcpp)(double xc, double yc, double *myxp)) {
      this->mapcpp = mapcpp;
    }
    //source updater
    void getSrcvalatIJ(int i, int j, double sitrij);
    //source updater
    void getConvalatIJ(int i, int j, int k, double citrij);
    //solution updater
    double getSolvalatIJ(int i, int j);

  private:
    //create laplacian (organize blocks)
    void laplace();
    //convert rhs source (only laplace currently)
    Eigen::VectorXd srcconv();
    //create blocks
    blkstor blkrow(int j);
    //coefficients
    double cij(double xc, double yc, int i, int j);
    double cimj(double xc, double yc, int i, int j);
    double cimjp(double xc, double yc, int i, int j);
    double cijp(double xc, double yc, int i, int j);
    double cipjp(double xc, double yc, int i, int j);
    double cipj(double xc, double yc, int i, int j);
    double cipjm(double xc, double yc, int i, int j);
    double cijm(double xc, double yc, int i, int j);
    double cimjm(double xc, double yc, int i, int j);
    //metric calculator
    void (*gfunc)(double *, double *);
    //diff length calculator
    void (*hfunc)(double *, double *);
    //mapc2p for special metric
    void (*mapcpp)(double, double, double *);
    //convert between gfunc output to eigen vector
    Eigen::Vector3d gij(double xc, double yc, int i, int j);
    //metric inverter
    Eigen::Vector3d ginv(Eigen::Vector3d gmat, double xc, double yc);
    //metric determinant
    double gdet(Eigen::Vector3d gvec);

    //block storage
    blkstor triblock;
    //simulation parameters
    int nx, ny, xbctypeu, xbctypel, ybctypeu, ybctypel;
    double xbcu, xbcl, ybcu, ybcl, xl, yl, xu, yu;
    //conductivity
    Eigen::MatrixXd condarr;
    //calculated parametrs
    double dxc, dyc;
    Eigen::MatrixXd xca;
    Eigen::MatrixXd yca;
    //solver type
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    //source array
    Eigen::VectorXd jsource;
    //solution array
    Eigen::VectorXd phi;
    //regular poisson or tensor Poisson
    bool sigma;
};

#endif
