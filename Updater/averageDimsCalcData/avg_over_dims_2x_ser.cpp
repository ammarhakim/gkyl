#include <avg_over_dims_mod_decl.h>

void avg_over_dims_2x_p1_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[2]*rNumCells; 

}

void avg_over_dims_2x_p1_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[1]*rNumCells; 

}


void avg_over_dims_2x_p2_ser_avgDirs1(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[2]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[5]*rNumCells; 

}

void avg_over_dims_2x_p2_ser_avgDirs2(double rNumCells, const double *fIn, double *fAvgOut) 
{
  // rNumCells: reciprocal of the number of cells over which we will average.
  // fIn: input field to average.
  // fAvgOut: output field average.

  fAvgOut[0] += 0.7071067811865475*fIn[0]*rNumCells; 
  fAvgOut[1] += 0.7071067811865475*fIn[1]*rNumCells; 
  fAvgOut[2] += 0.7071067811865475*fIn[4]*rNumCells; 

}


