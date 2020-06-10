// Gkyl ------------------------------------------------------------------------
//
// Maxwell equation object
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylCartField.h>
#include <GkylEquation.h>
#include <MaxwellTmplModDecl.h>
#include <GkylBasisTypes.h>

// std includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

extern "C" {

    typedef double (*Maxwell_volumeTerm_t)(const MaxwellEq_t *meq, const double*  w, const double*  dx,
      const double* q, double *out);
    typedef double (*Maxwell_surfTerm_t)(unsigned dir, const MaxwellEq_t *meq, const double*  wl, const double*  wr,
      const double*  dxl, const double*  dxr,
      const double tau,
      const double *ql, const double*  qr,
      double *outl, double *outr);
    
    typedef struct {
        // dims, basis info
        unsigned cdim, polyOrder, basisType;
        MaxwellEq_t mdata;
        double tau;
        
        // Maxwell-specific function pointers
        Maxwell_volumeTerm_t volumeTerm;
        Maxwell_surfTerm_t surfTerm;
    } GkylMaxwell_t;

    // Return a pointer to an equation object for Maxwell equations
    GkylEquation_t *new_MaxwellOnDevice(unsigned cdim, unsigned polyOrder, unsigned basisType,
      MaxwellEq_t *mdata, double tau);
}

