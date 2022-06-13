#ifndef CARTFIELD_AVERAGE_OVER_DIMS_MOD_DECL_H
#define CARTFIELD_AVERAGE_OVER_DIMS_MOD_DECL_H
extern "C" { 

void avg_over_dims_2x_p1_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_2x_p1_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 

void avg_over_dims_2x_p2_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_2x_p2_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 


void avg_over_dims_3x_p1_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p1_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p1_ser_avgDirs3(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p1_ser_avgDirs12(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p1_ser_avgDirs13(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p1_ser_avgDirs23(double rNumCells, const double *fIn, double *fAvgOut); 

void avg_over_dims_3x_p2_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p2_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p2_ser_avgDirs3(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p2_ser_avgDirs12(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p2_ser_avgDirs13(double rNumCells, const double *fIn, double *fAvgOut); 
void avg_over_dims_3x_p2_ser_avgDirs23(double rNumCells, const double *fIn, double *fAvgOut); 





} 

#endif 
