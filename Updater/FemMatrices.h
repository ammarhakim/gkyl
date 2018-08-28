#include <Eigen/Core> 
void getMassMatrix(Eigen::MatrixXd& mat, double *w, int ndim, int polyOrder);
void getMassMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder);
void getPerpNodToModMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder);
void getParNodToModMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder);
void getPerpStiffnessMatrix(Eigen::MatrixXd& mat, double *w, int ndim, int polyOrder, double dx, double dy);
void getPerpStiffnessMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder, double dx, double dy);
void getParStiffnessMatrix(Eigen::MatrixXd& mat, double *w, int ndim, int polyOrder, double dz);
void getParStiffnessMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder, double dz);
void getPerpStiffnessMatrix(Eigen::MatrixXd& mat, double *w, double *gxx, double *gxy, double *gyy, int ndim, int polyOrder, double dx, double dy);
void getModalMassMatrix(Eigen::MatrixXd& mat, double *w, int ndim, int polyOrder);
