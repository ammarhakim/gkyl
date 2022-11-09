-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the primitive moments, u and vtSq=sqrt(T/m), given
-- the distribution function and its moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase     = require "Updater.Base"
local Lin             = require "Lib.Linalg"
local Proto           = require "Lib.Proto"
local PrimMomentsDecl = require "Updater.primMomentsCalcData.PrimMomentsModDecl"
local xsys            = require "xsys"
local ffi             = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[
// Object type
typedef struct gkyl_prim_lbo_calc gkyl_prim_lbo_calc;

/**
 * Compute primitive moments of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param cbasis_rng Config-space basis functions
 * @param conf_rng Config-space range
 * @param moms Moments of distribution function (Zeroth, First, and Second)
 * @param boundary_corrections Momentum and Energy boundary corrections
 * @param prim_moms_out Output drift velocity and thermal speed squared.
 */
void gkyl_prim_lbo_calc_advance(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array *prim_moms_out);

void gkyl_prim_lbo_calc_advance_cu(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array* prim_moms_out);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* calc);

// "derived" class constructors
gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_with_fluid_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_with_fluid_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng);


// Object type
typedef struct gkyl_mom_calc_bcorr gkyl_mom_calc_bcorr;

/**
 * Compute boundary correction moments.
 *
 * @param bcorr Boundary correction updater object
 * @param phase_rng Phase space range on which to compute.
 * @param conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param out Output
 */
void gkyl_mom_calc_bcorr_advance(gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);

void gkyl_mom_calc_bcorr_advance_cu(gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param bcorr Updater to delete.
 */
void gkyl_mom_calc_bcorr_release(gkyl_mom_calc_bcorr* up);

// "derived" class constructors
gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary);

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary, double mass);

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary);

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_gyrokinetic_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary, double mass);
]]

-- function to check if operator option is correct
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" then
      return true
   end
   return false
end

-- Moments updater object.
local SelfPrimMoments = Proto(UpdaterBase)

function SelfPrimMoments:init(tbl)
   SelfPrimMoments.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SelfPrimMoments: Must provide grid object using 'onGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.SelfPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   self.confBasis = assert(
      tbl.confBasis, "Updater.SelfPrimMoments: Must provide the configuration basis object using 'confBasis'.")
   local confBasis = self.confBasis

   local operator = assert(
      tbl.operator, "Updater.SelfPrimMoments: Must specify the collision operator (VmLBO or GkLBO for now) using 'operator'.")

   -- Ensure sanity.
   local basisID = confBasis:id()  -- Basis name.
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
          "Polynomial orders of phase and conf basis must match.")
   assert((phaseBasis:id() == basisID) or
          ((phaseBasis:id()=="hybrid" or phaseBasis:id()=="gkhybrid") and basisID=="serendipity"),
          "Type of phase and conf basis must match.")

   self._isGkLBO    = false
   if isOperatorGood(operator) then
     if operator == "GkLBO" then self._isGkLBO = true end
   else
      assert(false, string.format(
                "SelfPrimMoments: Operator option must be 'VmLBO' or 'GkLBO'. Requested %s instead.", operator))
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   local vbounds = assert(tbl.vbounds, "SelfPrimMoments: must pass vbounds")
   local mass
   if self._isGkLBO then mass = assert(tbl.mass, "SelfPrimMoments: must pass species mass for GkLBO") end
   if GKYL_USE_GPU then
      if self._isGkLBO then
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_gyrokinetic_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds, mass),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_calc_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      else
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_vlasov_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_calc_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      end
   else
      if self._isGkLBO then
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds, mass),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_calc_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      else
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_vlasov_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_calc_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      end
   end
end

-- Advance method.
function SelfPrimMoments:_advance(tCurr, inFld, outFld)

   local moments, fIn = inFld[1], inFld[2]
   local boundaryCorrections, primMoms = outFld[1], outFld[2]

   -- Compute boundary corrections.
   ffiC.gkyl_mom_calc_bcorr_advance(self._zero_bcorr_calc, fIn:localRange(), primMoms:localRange(), fIn._zero, boundaryCorrections._zero)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_lbo_calc_advance(self._zero_prim_calc, self.confBasis._zero, primMoms:localRange(), moments._zero, boundaryCorrections._zero, primMoms._zero)

end

function SelfPrimMoments:_advanceOnDevice(tCurr, inFld, outFld)

   local moments, fIn = inFld[1], inFld[2]
   local boundaryCorrections, primMoms = outFld[1], outFld[2]

   -- Compute boundary corrections.
   ffiC.gkyl_mom_calc_bcorr_advance_cu(self._zero_bcorr_calc, fIn:localRange(), primMoms:localRange(), fIn._zeroDevice, boundaryCorrections._zeroDevice)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_lbo_calc_advance_cu(self._zero_prim_calc, self.confBasis._zero, primMoms:localRange(), moments._zeroDevice, boundaryCorrections._zeroDevice, primMoms._zeroDevice)

end
   

return SelfPrimMoments
