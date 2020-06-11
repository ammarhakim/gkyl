// Gkyl ------------------------------------------------------------------------
//
// Vlasov equation object
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylCartField.h>
#include <GkylEquation.h>
#include <VlasovTmplModDecl.h>
#include <GkylBasisTypes.h>

// std includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

extern "C" {

    typedef double (*Vlasov_volumeStreamTerm_t)(const double* __restrict__  w, const double* __restrict__  dxv, const double* __restrict__  f, double *out);
    typedef void (*Vlasov_surfSreamTerm_t)(unsigned dir, const double* __restrict__  wl, const double* __restrict__  wr,
      const double* __restrict__  dxvl, const double* __restrict__  dxvr, const double* __restrict__  fl, const double* __restrict__  fr,
      double *outl, double *outr);

    typedef double (*Vlasov_volumeTerm_t)(const double* __restrict__  w, const double* __restrict__  dxv,
      const double* __restrict__  E, const double* __restrict__  f, double *out);
    typedef double (*Vlasov_surfElcMagTerm_t)(unsigned dir, const double* __restrict__  wl, const double* __restrict__  wr,
      const double* __restrict__  dxvl, const double* __restrict__  dxvr,
      const double amax, const double* __restrict__  E, const
      double *fl, const double* __restrict__  fr,
      double *outl, double *outr);

    typedef struct {
        // dims, basis info
        unsigned cdim, vdim, polyOrder, basisType;
        // species parameters
        double qbym;
        bool hasForceTerm;
        // pointer to EM field
        GkylCartField_t *emField;

        // Vlasov-specific function pointers
        Vlasov_volumeStreamTerm_t volumeStreamTerm;
        Vlasov_surfSreamTerm_t surfStreamTerm;
        Vlasov_volumeTerm_t volumeTerm;
        Vlasov_surfElcMagTerm_t surfElcMagTerm;

    } GkylVlasov_t;

    // Return a pointer to an equation object for Vlasov equations
    GkylEquation_t *new_VlasovOnDevice(unsigned cdim, unsigned vdim, unsigned polyOrder, unsigned basisType,
      double qbym, bool hasForceTerm);
    // Set the aux fields
    void Vlasov_setAuxFields(GkylEquation_t *eqn, GkylCartField_t* em);
}
