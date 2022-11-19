-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create the Maxwellian distribution from the conserved
-- moments and project it on basis functions. Uses Gaussian
-- quadrature.
--
-- Note: There's an implementation below that has more calculations in Lua.
--       This is slower but it is kept because as of 11/03/2020 it is the only
--       way to do quadrature with >polyOrder+1 quad points. We could add support
--       for such quadratures in the C implementation by generating more kernels.
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
typedef struct gkyl_proj_maxwellian_on_basis gkyl_proj_maxwellian_on_basis;

/**
 * Create new updater to project Maxwellian on basis functions. Free
 * using gkyl_proj_maxwellian_on_basis_release method.
 *
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes (in 1D).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_proj_maxwellian_on_basis* gkyl_proj_maxwellian_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, bool use_gpu);

/**
 * Compute projection of Maxwellian on basis. This method takes
 * lab-frame moments to compute the projection of Maxwellian on basis
 * functions.
 *
 * @param mob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param fmax Output Maxwellian
 */
void gkyl_proj_maxwellian_on_basis_lab_mom(const gkyl_proj_maxwellian_on_basis *mob,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms, struct gkyl_array *fmax);

/**
 * Compute projection of Maxwellian on basis. This method takes
 * primitive (fluid-frame) moments to compute the projection of
 * Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param prim_moms (primitive moments udrift, vtsq=T/m)
 * @param fmax Output Maxwellian
 */
void gkyl_proj_maxwellian_on_basis_prim_mom(const gkyl_proj_maxwellian_on_basis *mob,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms, const struct gkyl_array *prim_moms, struct gkyl_array *fmax);

/**
 * Compute projection of a gyrokinetic Maxwellian on basis.
 * This method takes lab-frame moments to compute the projection
 * of Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param bmag Magnetic field magnitude.
 * @param jacob_tot Total jacobian (conf * guiding center jacobian).
 * @param mass Species mass.
 * @param fmax Output Maxwellian
 */
void gkyl_proj_gkmaxwellian_on_basis_lab_mom(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *bmag,
  const struct gkyl_array *jacob_tot, double mass, struct gkyl_array *fmax);

/**
 * Compute projection of a gyrokinetic Maxwellian on basis. This
 * method takes primitive (fluid-frame) moments to compute the
 * projection of Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param prim_moms (primitive moments upar, vtsq=T/m)
 * @param bmag Magnetic field magnitude.
 * @param jacob_tot Total jacobian (conf * guiding center jacobian).
 * @param mass Species mass.
 * @param fmax Output Maxwellian
 */
void gkyl_proj_gkmaxwellian_on_basis_prim_mom(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *prim_moms,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass,
  struct gkyl_array *fmax);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_maxwellian_on_basis_release(gkyl_proj_maxwellian_on_basis* mob);
]]


-- Inherit the base Updater from UpdaterBase updater object.
local MaxwellianOnBasis = Proto(UpdaterBase)

function MaxwellianOnBasis:init(tbl)
   MaxwellianOnBasis.super.init(self, tbl) -- setup base object

   local phaseGrid  = assert(tbl.onGrid,
                            "Updater.MaxwellianOnBasis: Must provide phase space grid object 'onGrid'")
   local confBasis  = assert(tbl.confBasis,
                            "Updater.MaxwellianOnBasis: Must provide configuration space basis object 'confBasis'")
   local phaseBasis = assert(tbl.phaseBasis,
                            "Updater.MaxwellianOnBasis: Must provide phase space basis object 'phaseBasis'")

   self.mass = tbl.mass  -- Mass needed to project a gyrokinetic Maxwellian.
   local usePrimMoms = xsys.pickBool(tbl.usePrimMoms, false)

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   -- Number of quadrature points in each direction.
   local numQuad1D = confBasis:polyOrder()+1
   self._zero = ffi.gc(ffiC.gkyl_proj_maxwellian_on_basis_new(phaseGrid._zero, confBasis._zero, 
                                                              phaseBasis._zero, numQuad1D, self._useGPU),
                       ffiC.gkyl_proj_maxwellian_on_basis_release)

   local isGK = false
   if self.mass then isGK = true end

   if isGK then
      self.advanceFunc = usePrimMoms and
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advance_prim_mom_gk"](self, tCurr, inFld, outFld)
         end or
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advance_lab_mom_gk"](self, tCurr, inFld, outFld)
         end
      self.advanceOnDeviceFunc = usePrimMoms and
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advanceOnDevice_prim_mom_gk"](self, tCurr, inFld, outFld)
         end or
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advanceOnDevice_lab_mom_gk"](self, tCurr, inFld, outFld)
         end
   else
      self.advanceFunc = usePrimMoms and 
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advance_prim_mom"](self, tCurr, inFld, outFld)
         end or
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advance_lab_mom"](self, tCurr, inFld, outFld)
         end
      self.advanceOnDeviceFunc = usePrimMoms and 
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advanceOnDevice_prim_mom"](self, tCurr, inFld, outFld)
         end or
         function(tCurr, inFld, outFld)
            MaxwellianOnBasis["_advanceOnDevice_lab_mom"](self, tCurr, inFld, outFld)
         end
   end

end

function MaxwellianOnBasis:_advance_lab_mom(tCurr, inFld, outFld)
   local momsIn = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local fOut   = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_maxwellian_on_basis_lab_mom(self._zero, phaseRange, confRange, momsIn._zero, fOut._zero)
end

function MaxwellianOnBasis:_advance_prim_mom(tCurr, inFld, outFld)
   local momsIn     = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local primMomsIn = assert(inFld[2], "MaxwellianOnBasis: Must specify primitive moments in 'inFld[2]'")
   local fOut       = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_maxwellian_on_basis_prim_mom(self._zero, phaseRange, confRange, momsIn._zero, primMomsIn._zero, fOut._zero)
end

function MaxwellianOnBasis:_advance_lab_mom_gk(tCurr, inFld, outFld)
   local momsIn   = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local bmag     = assert(inFld[2], "MaxwellianOnBasis: Must specify the magnetic field magnitude in 'inFld[2]'")
   local jacobTot = assert(inFld[3], "MaxwellianOnBasis: Must specify the total jacobian in 'inFld[2]'")
   local fOut     = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_gkmaxwellian_on_basis_lab_mom(self._zero, phaseRange, confRange, momsIn._zero,
                                                bmag._zero, jacobTot._zero, self.mass, fOut._zero)
end

function MaxwellianOnBasis:_advance_prim_mom_gk(tCurr, inFld, outFld)
   local momsIn     = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local primMomsIn = assert(inFld[2], "MaxwellianOnBasis: Must specify the primitive moments in 'inFld[2]'")
   local bmag       = assert(inFld[3], "MaxwellianOnBasis: Must specify the magnetic field magnitude in 'inFld[3]'")
   local jacobTot   = assert(inFld[4], "MaxwellianOnBasis: Must specify the total jacobian in 'inFld[4]'")
   local fOut       = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_gkmaxwellian_on_basis_prim_mom(self._zero, phaseRange, confRange, momsIn._zero, primMomsIn._zero,
                                                 bmag._zero, jacobTot._zero, self.mass, fOut._zero)
end

function MaxwellianOnBasis:_advanceOnDevice_lab_mom(tCurr, inFld, outFld)
   local momsIn = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local fOut   = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_maxwellian_on_basis_lab_mom(self._zero, phaseRange, confRange, momsIn._zeroDevice, fOut._zeroDevice)
end

function MaxwellianOnBasis:_advanceOnDevice_prim_mom(tCurr, inFld, outFld)
   local momsIn     = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local primMomsIn = assert(inFld[2], "MaxwellianOnBasis: Must specify primitive moments in 'inFld[2]'")
   local fOut       = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_maxwellian_on_basis_prim_mom(self._zero, phaseRange, confRange, momsIn._zeroDevice,
                                               primMomsIn._zeroDevice, fOut._zeroDevice)
end

function MaxwellianOnBasis:_advanceOnDevice_lab_mom_gk(tCurr, inFld, outFld)
   local momsIn   = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local bmag     = assert(inFld[2], "MaxwellianOnBasis: Must specify the magnetic field magnitude in 'inFld[2]'")
   local jacobTot = assert(inFld[3], "MaxwellianOnBasis: Must specify the total jacobian in 'inFld[2]'")
   local fOut     = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_gkmaxwellian_on_basis_lab_mom(self._zero, phaseRange, confRange, momsIn._zeroDevice,
                                                bmag._zeroDevice, jacobTot._zeroDevice, self.mass, fOut._zeroDevice)
end

function MaxwellianOnBasis:_advanceOnDevice_prim_mom_gk(tCurr, inFld, outFld)
   local momsIn     = assert(inFld[1], "MaxwellianOnBasis: Must specify the velocity moments in 'inFld[1]'")
   local primMomsIn = assert(inFld[2], "MaxwellianOnBasis: Must specify the primitive moments in 'inFld[2]'")
   local bmag       = assert(inFld[3], "MaxwellianOnBasis: Must specify the magnetic field magnitude in 'inFld[3]'")
   local jacobTot   = assert(inFld[4], "MaxwellianOnBasis: Must specify the total jacobian in 'inFld[4]'")
   local fOut       = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field in 'outFld[1]'")

   local confRange  = self.onGhosts and momsIn:localExtRange() or momsIn:localRange()
   local phaseRange = self.onGhosts and fOut:localExtRange() or fOut:localRange()

   ffiC.gkyl_proj_gkmaxwellian_on_basis_prim_mom(self._zero, phaseRange, confRange, momsIn._zeroDevice, primMomsIn._zeroDevice,
                                                 bmag._zeroDevice, jacobTot._zeroDevice, self.mass, fOut._zeroDevice)
end

function MaxwellianOnBasis:_advance(tCurr, inFld, outFld)
   self.advanceFunc(tCurr, inFld, outFld)
end

function MaxwellianOnBasis:_advanceOnDevice(tCurr, inFld, outFld)
   self.advanceOnDeviceFunc(tCurr, inFld, outFld)
end

return MaxwellianOnBasis
