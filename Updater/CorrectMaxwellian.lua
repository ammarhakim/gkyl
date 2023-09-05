-- Gkyl ------------------------------------------------------------------------
--
-- Updater to conserve the moments of the Maxwellian with the iteration fix.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto = require "Proto"
local UpdaterBase = require "Updater.Base"
local xsys = require "xsys"
local ffi = require "ffi"

local ffiC = ffi.C
ffi.cdef[[
// Object type
typedef struct gkyl_correct_maxwellian_gyrokinetic gkyl_correct_maxwellian_gyrokinetic;

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_local Local configuration space range
 * @param conf_local_ext Local extended configuration space range
 * @param bmag Magnetic field
 * @param jacob_tot Jacobian to project the maxwellian distribution
 * @param mass Mass of the species
 * @param use_gpu Bool to determine if on GPU
 */
gkyl_correct_maxwellian_gyrokinetic *gkyl_correct_maxwellian_gyrokinetic_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass, bool use_gpu);

/**
 * Fix the Maxwellian so that it's moments match desired moments.
 *
 * @param cmax Maxwellian-fix updater
 * @param fM Distribution function to fix (modified in-place)
 * @param moms_in Input moments
 * @param err_max Tolerance of error in M1 and M2
 * @param iter_max Maximum number of iteration
 * @param conf_local Local configuration space range
 * @param conf_local_ext Local extended configuration space range
 * @param phase_local Local phase-space range
 */
void gkyl_correct_maxwellian_gyrokinetic_fix(gkyl_correct_maxwellian_gyrokinetic *cmax,
  struct gkyl_array *fM, const struct gkyl_array *moms_in, double err_max, int iter_max,
  const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext,
  const struct gkyl_range *phase_local);

/**
 * Delete updater.
 *
 * @param cmax Updater to delete.
 */
void gkyl_correct_maxwellian_gyrokinetic_release(gkyl_correct_maxwellian_gyrokinetic* cmax);
]]

-- Inherit the base Updater from UpdaterBase updater object.
local CorrectMaxwellian = Proto(UpdaterBase)

function CorrectMaxwellian:init(tbl)
   CorrectMaxwellian.super.init(self, tbl)  -- set up base object

   local phaseGrid = assert(tbl.onGrid, "Updater.CorrectMaxwellian: Must provide phase space grid object 'onGrid'")
   local confGrid = assert(tbl.confGrid, "Updater.CorrectMaxwellian: Must provide configuration space grid object 'confGrid'")
   local phaseBasis = assert(tbl.phaseBasis, "Updater.CorrectMaxwellian: Must provide phase space basis object 'phaseBasis'")
   local confBasis = assert(tbl.confBasis, "Updater.CorrectMaxwellian: Must provide configuration space basis object 'confBasis'")
   local bmag = assert(tbl.bmag, "Updater.CorrectMaxwellian: Must provide gkyl_array object 'bmag'")
   local jacobTot = assert(tbl.jacobTot, "Updater.CorrectMaxwellian: Must provide gkyl_array object 'jacobTot'")
   
   --self.bmag = tbl.bmag
   --self.jacobTot = tbl.jacobTot
   self.mass = tbl.mass
   self.iter_max = tbl.iter_max
   self.err_max = tbl.err_max
   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)
 
   self.confRange = bmag:localRange() 
   self.confRangeExt = bmag:localExtRange()

   self._zero = ffi.gc(ffiC.gkyl_correct_maxwellian_gyrokinetic_new(phaseGrid._zero, confBasis._zero, phaseBasis._zero, self.confRange, self.confRangeExt, bmag._zero, jacobTot._zero, self.mass, self._useGPU), ffiC.gkyl_correct_maxwellian_gyrokinetic_release)
 
end

function CorrectMaxwellian:_advance(tCurr, inFld, outFld)
   local fM = assert(inFld[1], "CorrectMaxwellian.advance: Must specify an input maxwellian in 'inFld[1]'")
   local moms = assert(inFld[2], "CorrectMaxwellian.advance: Must specify the target velocity moments in 'inFld[2]'")
   local fOut =  assert(outFld[1], "CorrectMaxwellian.advance: Must specify an output field in 'outFld[1]'")

   local phaseRange = fM:localRange() 

   print("CorrectMaxwellian is called")

   ffiC.gkyl_correct_maxwellian_gyrokinetic_fix(self._zero, fM._zero, moms._zero, self.err_max, self.iter_max, self.confRange, self.confRangeExt, phaseRange)
   fOut:copy(fM)
end

--[[function CorrectMaxwellian:_advanceOnDevice(tCurr, inFld, outFld)
end--]]

return CorrectMaxwellian

