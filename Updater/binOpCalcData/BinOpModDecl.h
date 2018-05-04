#ifndef BIN_OP_MOD_DECL_H 
#define BIN_OP_MOD_DECL_H 
 
// Eigen include statements. 
#include <Eigen/Dense> 
 
extern "C" { 
void BinOpMultiply1xSer_P1(const double *A, const double *B, const short int vDim, double *out); 

void BinOpMultiply1xSer_P2(const double *A, const double *B, const short int vDim, double *out); 


void BinOpDivide1xSer_P1(const double *A, const double *B, const short int vDim, double *out); 

void BinOpDivide1xSer_P2(const double *A, const double *B, const short int vDim, double *out); 


} 
#endif 
