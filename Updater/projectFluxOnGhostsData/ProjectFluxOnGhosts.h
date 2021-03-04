#ifndef PROJECTFLUX_MOD_DECL_H 
#define PROJECTFLUX_MOD_DECL_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
void ProjectFluxOnGhosts1x1vDir1Ser_P1(const double zVal, const double *fIn, double *fHat); 

void ProjectFluxOnGhosts1x2vDir1Ser_P1(const double zVal, const double *fIn, double *fHat); 

void ProjectFluxOnGhosts1x3vDir1Ser_P1(const double zVal, const double *fIn, double *fHat); 


void ProjectFluxOnGhosts2x2vDir2Ser_P1(const double zVal, const double *fIn, double *fHat); 

void ProjectFluxOnGhosts2x3vDir2Ser_P1(const double zVal, const double *fIn, double *fHat); 


void ProjectFluxOnGhosts3x3vDir3Ser_P1(const double zVal, const double *fIn, double *fHat); 


} 
#endif 
