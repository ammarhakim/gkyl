#include <avg_over_dims_mod_decl.h>

void avg_over_dims_3x_p1_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[3]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[3]*rNumCells; 
  fAvgOut[3] += 0.7071067811865475*fIn[0]*rNumCells; 

}

void avg_over_dims_3x_p1_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[1]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[3]*rNumCells; 
  fAvgOut[3] += 0.7071067811865475*fIn[5]*rNumCells; 

}

void avg_over_dims_3x_p1_ser_avgDirs3(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[1]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[2]*rNumCells; 
  fAvgOut[3] += 0.7071067811865475*fIn[4]*rNumCells; 

}

void avg_over_dims_3x_p1_ser_avgDirs12(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.5*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.5*fIn[3]*rNumCells; 

}

void avg_over_dims_3x_p1_ser_avgDirs13(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.5*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.5*fIn[2]*rNumCells; 

}

void avg_over_dims_3x_p1_ser_avgDirs23(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.5*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.5*fIn[1]*rNumCells; 

}


void avg_over_dims_3x_p2_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[3]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[3]*rNumCells; 
  fAvgOut[3] += 0.6324555320336759*fIn[9]*rNumCells+0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[4] += 0.7071067811865475*fIn[9]*rNumCells; 
  fAvgOut[5] += 0.7071067811865475*fIn[9]*rNumCells; 
  fAvgOut[6] += 0.632455532033676*fIn[3]*rNumCells; 
  fAvgOut[7] += 0.632455532033676*fIn[3]*rNumCells; 

}

void avg_over_dims_3x_p2_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[1]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[3]*rNumCells; 
  fAvgOut[3] += 0.7071067811865475*fIn[5]*rNumCells; 
  fAvgOut[4] += 0.7071067811865475*fIn[7]*rNumCells; 
  fAvgOut[5] += 0.7071067811865475*fIn[9]*rNumCells; 
  fAvgOut[6] += 0.7071067811865475*fIn[13]*rNumCells; 
  fAvgOut[7] += 0.7071067811865475*fIn[15]*rNumCells; 

}

void avg_over_dims_3x_p2_ser_avgDirs3(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[1]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[2]*rNumCells; 
  fAvgOut[3] += 0.7071067811865475*fIn[4]*rNumCells; 
  fAvgOut[4] += 0.7071067811865475*fIn[7]*rNumCells; 
  fAvgOut[5] += 0.7071067811865475*fIn[8]*rNumCells; 
  fAvgOut[6] += 0.7071067811865475*fIn[11]*rNumCells; 
  fAvgOut[7] += 0.7071067811865475*fIn[12]*rNumCells; 

}

void avg_over_dims_3x_p2_ser_avgDirs12(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.5*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.5*fIn[3]*rNumCells; 
  fAvgOut[2] += 0.5*fIn[9]*rNumCells; 

}

void avg_over_dims_3x_p2_ser_avgDirs13(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.5*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.5*fIn[2]*rNumCells; 
  fAvgOut[2] += 0.5*fIn[8]*rNumCells; 

}

void avg_over_dims_3x_p2_ser_avgDirs23(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.5*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.5*fIn[1]*rNumCells; 
  fAvgOut[2] += 0.5*fIn[7]*rNumCells; 

}


