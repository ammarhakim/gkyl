#ifndef CARTFIELD_AVERAGE_OVER_DIMS_MOD_DECL_H
#define CARTFIELD_AVERAGE_OVER_DIMS_MOD_DECL_H
extern "C" { 

double avg_over_dims_2x_p1_Ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_2x_p1_Ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 

double avg_over_dims_2x_p2_Ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_2x_p2_Ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 


double avg_over_dims_3x_p1_Ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p1_Ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p1_Ser_avgDirs3(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p1_Ser_avgDirs12(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p1_Ser_avgDirs13(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p1_Ser_avgDirs23(double rNumCells, const double *fIn, double *fAvgOut); 

double avg_over_dims_3x_p2_Ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p2_Ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p2_Ser_avgDirs3(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p2_Ser_avgDirs12(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p2_Ser_avgDirs13(double rNumCells, const double *fIn, double *fAvgOut); 
double avg_over_dims_3x_p2_Ser_avgDirs23(double rNumCells, const double *fIn, double *fAvgOut); 





} 

#endif 
