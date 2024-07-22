-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate special relativistic variables for use in 
-- special relativistic Vlasov simulations.
-- Currently supporting the following:
--   1) p/gamma, gamma, & 1/gamma initialization, momentum-space grid objects where gamma = sqrt(1 + p^2/c^2).
--   2) GammaV2 = 1/(1 - V^2/c^2), square of the Lorentz boost factor for a given bulk velocity, V
--   3) GammaV = 1/sqrt(1 - V^2/c^2), Lorentz boost factor for a given bulk velocity, V.
--   4) GammaV_inv = sqrt(1 - V^2/c^2), inverse of the Lorentz boost factor for a given bulk velocity, V.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local xsys        = require "xsys"
local ffi         = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[
// Object type
typedef struct gkyl_dg_calc_sr_vars gkyl_dg_calc_sr_vars;

/**
 * Create new updater to compute relativistic variables needed in 
 * updates and used for diagnostics. Methods compute:
 * "p_vars" : the particle Lorentz boost factor gamma = sqrt(1 + p^2) and its inverse
 * n : the rest-frame plasma density n = GammaV_inv*M0 where GammaV_inv = sqrt(1 - |V_drift|^2)
 *     and V_drift is the bulk velocity computed via weak division V_drift = M1i/M0
 * u_i : the bulk four-velocity (GammaV, GammaV*V_drift) computed using weak division from M0 and M1i
 * p = n*T : the rest-frame pressure. The rest-frame pressure is computed as a velocity moment 
 *           with the weight = gamma*GammaV^2 - 2*GammaV*(v . p) + 1/gamma*((v . p)^2 - 1)
 *           where v is the spatial component of the bulk four-velocity: GammaV*V_drift, 
 *           GammaV is the bulk Lorentz boost factor: sqrt(1 + v^2), 
 *           p is the spatial component of the particle four-velocity, 
 *           and gamma = sqrt(1 + p^2) is the particle Lorentz boost factor.
 * 
 * @param phase_grid Phase-space grid (for getting cell spacing and cell center in pressure calculation) 
 * @param vel_grid   Momentum (four-velocity)-space grid 
 *                   (for getting cell spacing and cell center in gamma and 1/gamma calculation) 
 * @param conf_basis Configuration-space basis functions
 * @param vel_basis  Momentum (four-velocity)-space basis functions
 * @param mem_range  Configuration-space range that sets the size of the bin_op memory
 *                   for computing V_drift. Note range is stored so updater loops 
 *                   over consistent range solving linear systems since memory is pre-allocated.
 * @param vel_range  Momentum (four-velocity)-space
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_sr_vars* 
gkyl_dg_calc_sr_vars_new(const struct gkyl_rect_grid *phase_grid, const struct gkyl_rect_grid *vel_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *vel_basis, 
  const struct gkyl_range *mem_range, const struct gkyl_range *vel_range, bool use_gpu);

/**
 * Compute the momentum grid variables for special relativistic simulations
 * Uses special kernels which convert between a Gauss-Lobatto nodal basis and
 * our modal basis to insure continuity of the momentum (four-velocity)-grid variables.
 *
 * @param up        Updater for computing sr variables 
 * @param gamma     Output array of particle Lorentz boost factor, gamma = sqrt(1 + p^2) 
 * @param gamma_inv Output array of inverse particle Lorentz boost factor, 1/gamma = 1/sqrt(1 + p^2) 
 */
void gkyl_calc_sr_vars_init_p_vars(struct gkyl_dg_calc_sr_vars *up, 
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv);

/**
 * Delete pointer to updater to compute sr variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_sr_vars_release(struct gkyl_dg_calc_sr_vars *up);  
]]

-- Moments updater object.
local CalcSRVars = Proto(UpdaterBase)

function CalcSRVars:init(tbl)
   CalcSRVars.super.init(self, tbl) -- Setup base object.

   -- Read data from input file.
   self._onGrid = assert(
      tbl.onGrid, "Updater.CalcSRVars: Must provide phase space grid object using 'onGrid'")
   self._velGrid = assert(
      tbl.velGrid, "Updater.CalcSRVars: Must provide velocity grid object using 'velGrid'")

   self._confBasis = assert(
      tbl.confBasis, "Updater.CalcSRVars: Must specify configuration space basis with 'confBasis'.")
   self._velBasis = assert(
      tbl.velBasis, "Updater.CalcSRVars: Must specify velocity space basis with 'velBasis'.")

   self._confRange = assert(
     tbl.confRange, "Updater.CalcSRVars: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Updater.CalcSRVars: confRange must be a sub-range") 
   self._velRange = assert(
     tbl.velRange, "Updater.CalcSRVars: Must specify velocity-space range using 'velRange'")
   assert(self._velRange:isSubRange()==1, "Updater.CalcSRVars: velRange must be a sub-range") 

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   self._zero = ffi.gc(
      ffiC.gkyl_dg_calc_sr_vars_new(self._onGrid._zero, self._velGrid._zero, 
        self._confBasis._zero, self._velBasis._zero, 
        self._confRange, self._velRange, self._useGPU),
      ffiC.gkyl_dg_calc_sr_vars_release
   )

   return self
end

-- Initialize SR variables (gamma, gamma_inv) advance method.
function CalcSRVars:_advance(tCurr, inFld, outFld)

   local gamma = outFld[1]
   local gamma_inv = outFld[2]

   ffiC.gkyl_calc_sr_vars_init_p_vars(self._zero, gamma._zero, gamma_inv._zero)
end

return CalcSRVars
