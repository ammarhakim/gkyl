-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function.
--
-- For LBO collisions this updater also computes the boundary corrections and,
-- if using a piecewise polynomial basis, the star moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local xsys        = require "xsys"
local ffi         = require "ffi"
local lume        = require "Lib.lume"

local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[ 
// Forward declare for use in function pointers
struct gkyl_mom_type;

/**
 * Function pointer type to compute the needed moment.
 */
typedef void (*momf_t)(const struct gkyl_mom_type *momt,
  const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param);

struct gkyl_mom_type {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  int num_mom; // number of components in moment
  momf_t kernel; // moment calculation kernel
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flag;
  struct gkyl_mom_type *on_dev; // pointer to itself or device data
};

/**
 * Create new Vlasov moment type object. Valid 'mom' strings are "M0",
 * "M1i", "M2", "M2ij", "M3i", "M3ijk", "FiveMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);

struct gkyl_mom_type* gkyl_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);

/**
 * Create new Gyrokinetic moment type object. Valid 'mom' strings are "GkM0",
 * "GkM1", "GkM2", "GkM2par", "GkM2perp", "GkM3par", "GkM3perp", "ThreeMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_range Configuration-space range
 * @param mass Mass of species
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, const char *mom);

/**
 * Create new Gyrokinetic moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, const char *mom);

/**
 * Set magnitude of the magnetic field, bmag, needed in computing moments.
 * 
 * @param momt Moment type pointer
 * @param bmag Pointer to magnitude of magnetic field
 */
void gkyl_gyrokinetic_set_bmag(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag);

struct gkyl_mom_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_mom_type *momt;

  uint32_t flags;
  struct gkyl_mom_calc *on_dev; // pointer to itself or device data
};

// Object type
typedef struct gkyl_mom_calc gkyl_mom_calc;

/**
 * Create new updater to compute moments of distribution
 * function. Free using gkyl_mom_calc_new_release.
 *
 * @param grid Grid object
 * @param momt Pointer to moment type object
 * @return New updater pointer.
 */
gkyl_mom_calc* gkyl_mom_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Create new updater to compute moments of distribution function on
 * NV-GPU. See new() method for documentation.
 */
gkyl_mom_calc* gkyl_mom_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Compute moment of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Moment calculator updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 * @param mout Output moment array
 */
void gkyl_mom_calc_advance(const gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout);

void gkyl_mom_calc_advance_cu(const gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout);
]]

-- Moments updater object.
local DistFuncMomentCalc = Proto(UpdaterBase)

-- Valid moment names for Vlasov and GK equations.
local goodMomNames = {
   "M0", "M1i", "M2ij", "M2", "M3i", "FiveMoments", "FiveMomentsLBO",
}
local goodGkMomNames = {
   "GkM0", "GkM1", "GkM1proj", "GkM2par", "GkM2perp", "GkM2", "GkM3par", "GkM3perp",
   "GkThreeMoments", "GkThreeMomentsLBO"
}
local goodPartialMomNames = {
   -- Partial velocity moments, integrating over the region
   -- where one of the velocities if positive or negative.
   "M0Pvx",  "M0Pvy",  "M0Pvz", "M0Nvx",  "M0Nvy",  "M0Nvz",
   "M1iPvx", "M1iPvy", "M1iPvz","M1iNvx", "M1iNvy", "M1iNvz",
   "M2Pvx",  "M2Pvy",  "M2Pvz", "M2Nvx",  "M2Nvy",  "M2Nvz",
   "M3iPvx", "M3iPvy", "M3iPvz","M3iNvx", "M3iNvy", "M3iNvz"
}

function DistFuncMomentCalc:isMomentNameGood(nm)
   if lume.find(goodMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:isGkMomentNameGood(nm)
   if lume.find(goodGkMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:isPartialMomentNameGood(nm)
   if lume.find(goodPartialMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:init(tbl)
   DistFuncMomentCalc.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncMomentCalc: Must provide grid object using 'onGrid'")
   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   local advArgs = tbl.advanceArgs  -- Sample arguments for advance method.

   local confBasisID, polyOrder = confBasis:id(), confBasis:polyOrder()

   -- Dimension of spaces.
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()

   -- Ensure sanity.
   assert(polyOrder == phaseBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert((confBasisID == phaseBasis:id()) or
          ((phaseBasis:id()=="hybrid" or phaseBasis:id()=="gkhybrid") and confBasisID=="serendipity"),
	  "Type of phase-space and config-space basis must match")

   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   local isPartialMom = false
   if mom == "FiveMoments" or mom == "GkThreeMoments" or
      mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then

      if mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
         -- Rename this variable to call the right kernel below.
         mom = mom=="FiveMomentsLBO" and "FiveMoments" or "GkThreeMoments"
      end
   elseif self:isPartialMomentNameGood(mom) then
      isPartialMom = true
      assert(false, "Updater.DistFuncMomentCalc: partial moments are not implemented.")
   end

   -- Function to compute specified moment.
   local isGk = false
   if self:isMomentNameGood(mom) then
   elseif self:isGkMomentNameGood(mom) then
      isGk = true
      assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table 
                           containing the species mass and the background magnetic field
                           to calculate a Gk moment]])
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i, FiveMoments, or FiveMomentsLBO")
   end

   local mass, bmag
   if tbl.gkfacs then
      mass = tbl.gkfacs[1]
      bmag = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, true)

   if GKYL_USE_GPU then
      if isGk then
         self._zero_mom_type = ffiC.gkyl_mom_gyrokinetic_cu_dev_new(confBasis._zero, phaseBasis._zero, bmag:localRange(), mass, string.sub(mom,3))
         ffiC.gkyl_gyrokinetic_set_bmag(self._zero_mom_type, bmag._zeroDevice)
      else
         self._zero_mom_type = ffiC.gkyl_mom_vlasov_cu_dev_new(confBasis._zero, phaseBasis._zero, mom)
      end
      self._zero = ffiC.gkyl_mom_calc_cu_dev_new(self._onGrid._zero, self._zero_mom_type)
   else
      if isGk then
         self._zero_mom_type = ffiC.gkyl_mom_gyrokinetic_new(confBasis._zero, phaseBasis._zero, bmag:localRange(), mass, string.sub(mom,3))
         ffiC.gkyl_gyrokinetic_set_bmag(self._zero_mom_type, bmag._zero)
      else
         self._zero_mom_type = ffiC.gkyl_mom_vlasov_new(confBasis._zero, phaseBasis._zero, mom)
      end
      self._zero = ffiC.gkyl_mom_calc_new(self._onGrid._zero, self._zero_mom_type)
   end
end

function DistFuncMomentCalc:_advance(tCurr, inFld, outFld)
   local distf, mout = inFld[1], outFld[1]
   local phaseRange = distf:localRange()
   -- need to use localExtRange for confRange when onGhosts=true
   -- note that phaseRange does not need to use localExtRange because
   -- phaseRange is only used to get the velocity-space sub-range,
   -- which should not include ghosts
   local confRange = self.onGhosts and mout:localExtRange() or mout:localRange()
   ffiC.gkyl_mom_calc_advance(self._zero, phaseRange, confRange, distf._zero, mout._zero)
end

function DistFuncMomentCalc:_advanceOnDevice(tCurr, inFld, outFld)
   local distf, mout = inFld[1], outFld[1]
   local phaseRange = distf:localRange()
   local confRange
   if self.onGhosts then
      -- on GPU, we do need to set up phaseRange with ghosts in conf 
      -- space but not in v-space as a subrange of phaseExtRange
      local phaseExtRange = distf:localExtRange()
      local lower = phaseExtRange:lowerAsVec()
      local upper = phaseExtRange:upperAsVec()
      for d=self._cDim+1, self._pDim do
         lower[d] = lower[d] + distf:lowerGhost()
         upper[d] = upper[d] - distf:upperGhost()
      end
      phaseRange = phaseExtRange:subRange(lower, upper)
      confRange = mout:localExtRange()
   else
      confRange = mout:localRange()
   end
   mout:clear(0.0)
   ffiC.gkyl_mom_calc_advance_cu(self._zero, phaseRange, confRange, distf._zeroDevice, mout._zeroDevice)
end

return DistFuncMomentCalc
