-- Gkyl ------------------------------------------------------------------------
--
-- Project the square root (to a power) of a CartField (or its reciprocal) onto
-- the basis using quadrature. So can compute (sqrt(f))^q where q can be negative.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto       = require "Proto"
local UpdaterBase = require "Updater.Base"
local xsys        = require "xsys"
local ffi         = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_proj_powsqrt_on_basis gkyl_proj_powsqrt_on_basis;

/**
 * Create new updater to project pow( sqrt(f), e ) onto the basis via quadrature.
 *
 * @param basis Basis object (configuration space).
 * @param num_quad Number of quadrature nodes.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_proj_powsqrt_on_basis* gkyl_proj_powsqrt_on_basis_new(
  const struct gkyl_basis *basis, int num_quad, bool use_gpu);

/**
 * Compute pow( sqrt(fIn), expIn) via quadrature.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space range
 * @param expIn Exponent.
 * @param fIn Input scalar field.
 * @param fOut Ouput scalar field.
 */
void gkyl_proj_powsqrt_on_basis_advance(const gkyl_proj_powsqrt_on_basis *up,
  const struct gkyl_range *range, double expIn, const struct gkyl_array *fIn,
  struct gkyl_array *fOut);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_powsqrt_on_basis_release(gkyl_proj_powsqrt_on_basis* up);
]]

-- Inherit the base Updater from UpdaterBase updater object.
local SqrtOnBasis = Proto(UpdaterBase)

function SqrtOnBasis:init(tbl)
   SqrtOnBasis.super.init(self, tbl) -- setup base object

   self.basis = assert(tbl.basis, "Updater.SqrtOnBasis: Must provide configuration space basis object 'basis'.")

   -- Number of quadrature points in each direction
   local numQuad1D = self.basis:polyOrder() + 1

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   self._zero = ffi.gc(ffiC.gkyl_proj_powsqrt_on_basis_new(self.basis._zero, numQuad1D, self._useGPU),
                       ffiC.gkyl_proj_powsqrt_on_basis_release) 
end

function SqrtOnBasis:_advance(tCurr, inFld, outFld)
   -- Get the inputs and outputs.
   local fIn  = assert(inFld[1], "SqrtOnBasis: Must specify an input field 'inFld[1]'.")
   local fOut = assert(outFld[1], "SqrtOnBasis: Must specify an output field 'outFld[1]'.")

   local exponent = inFld[2] or 1.0

   local range = self.onGhosts and fIn:localExtRange() or fIn:localRange()

   ffiC.gkyl_proj_powsqrt_on_basis_advance(self._zero, range, exponent, fIn._zero, fOut._zero)
end

function SqrtOnBasis:_advanceOnDevice(tCurr, inFld, outFld)
   -- Get the inputs and outputs.
   local fIn  = assert(inFld[1], "SqrtOnBasis: Must specify an input field 'inFld[1]'.")
   local fOut = assert(outFld[1], "SqrtOnBasis: Must specify an output field 'outFld[1]'.")

   local exponent = inFld[2] or 1.0

   local range = self.onGhosts and fIn:localExtRange() or fIn:localRange()

   ffiC.gkyl_proj_powsqrt_on_basis_advance(self._zero, range, exponent, fIn._zeroDevice, fOut._zeroDevice)
end

return SqrtOnBasis
