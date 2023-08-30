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
/**
 * Compute the momentum grid variables for special relativistic simulations
 * Uses project_on_basis with Gauss-Lobatto nodes to insure continuity of resulting modal projection
 *
 * @param vgrid Momentum-space grid
 * @param vbasis Momentum-space basis
 * @param vrange Momentum-space range
 * @param p_over_gamma Output array of relativistic velocity, v = p/(gamma) = p/sqrt(1 + p^2)
 * @param gamma Output array of particle Lorentz boost factor, gamma = sqrt(1 + p^2) 
 * @param gamma_inv Output array of inverse particle Lorentz boost factor, 1/gamma = 1/sqrt(1 + p^2) 
 */
void gkyl_calc_sr_vars_init_p_vars(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_basis *vbasis, const struct gkyl_range *vrange,
  struct gkyl_array* p_over_gamma, struct gkyl_array* gamma, struct gkyl_array* gamma_inv);

/**
 * Compute the square of the Lorentz boost factor for a given bulk velocity, V.
 * GammaV2 = 1/(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute (1 - V^2/c^2) using basis_exp_sq (see gkyl_basis_*_exp_sq.h in kernels/basis/)
 * 2. Compute 1/(1 - V^2/c^2) using basis_inv (see gkyl_basis_*_inv.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param Gamma2V Output array of the square of the Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma2(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2);

/**
 * Compute the Lorentz boost factor for a given bulk velocity, V.
 * GammaV = 1/sqrt(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute 1/(1 - V^2/c^2) using basis_exp_sq and basis_inv
 *    (see gkyl_basis_*_exp_sq.h and gkyl_basis_*_inv.h in kernels/basis/)
 * 2. Project onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param Gamma Output array of Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV);

/**
 * Compute the inverse of the Lorentz boost factor for a given bulk velocity, V.
 * GammaV_inv = sqrt(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute GammaV2_inv = 1 - V^2/c^2 using basis_exp_sq 
 *    (see gkyl_basis_*_exp_sq.h in kernels/basis/)
 * 2. Project onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param Gamma Output array of Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma_inv(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv);
]]

-- Function to check if moment name is correct.
local function isOpNameGood(nm)
   if nm == "init" or nm == "Gamma2" or nm=="Gamma" or nm == "Gamma_inv" then
      return true
   end
   return false
end

-- Moments updater object.
local CalcSRVars = Proto(UpdaterBase)

function CalcSRVars:init(tbl)
   CalcSRVars.super.init(self, tbl) -- Setup base object.

   -- Read data from input file.
   self._velGrid = assert(
      tbl.velGrid, "Updater.CalcSRVars: Must provide velocity grid object using 'velGrid'")

   self._confBasis = assert(
      tbl.confBasis, "Updater.CalcSRVars: Must specify configuration space basis with 'confBasis'.")
   self._velBasis = assert(
      tbl.velBasis, "Updater.CalcSRVars: Must specify velocity space basis with 'velBasis'.")
   self._phaseBasis = assert(
      tbl.phaseBasis, "Updater.CalcSRVars: Must specify phase space basis with 'phaseBasis'.")

   local op = assert(
      tbl.operation, "Updater.CalcSRVars: Must provide an operation using 'operation'.")

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)

   assert(isOpNameGood(op), string.format(
          "CalcSRVars: Operation must be one of init, Gamma2, Gamma, or Gamma_inv. Requested %s instead.", op))

   -- Set the advance method function.
   if op == "init" then
      self.advanceFunc = function(tCurr, inFlds, outFlds)
         CalcSRVars['_advanceInitSRVars'](self, tCurr, inFlds, outFlds) end
      self.advanceOnDeviceFunc = nil -- Involves proj_on_basis which is not on device
   elseif op == "Gamma2" then
      self.advanceFunc = function(tCurr, inFlds, outFlds)
         CalcSRVars['_advanceGamma2'](self, tCurr, inFlds, outFlds) end
      self.advanceOnDeviceFunc = nil -- No GPU implementation yet
   elseif op == "Gamma" then
      self.advanceFunc = function(tCurr, inFlds, outFlds)
         CalcSRVars['_advanceGamma'](self, tCurr, inFlds, outFlds) end
      self.advanceOnDeviceFunc = nil -- No GPU implementation yet
   elseif op == "Gamma_inv" then
      self.advanceFunc = function(tCurr, inFlds, outFlds)
         CalcSRVars['_advanceGammaInv'](self, tCurr, inFlds, outFlds) end
      self.advanceOnDeviceFunc = nil -- No GPU implementation yet
   end
end

-- Initialize SR variables (p/gamma, gamma, gamma_inv) advance method.
function CalcSRVars:_advanceInitSRVars(tCurr, inFld, outFld)

   local p_over_gamma = outFld[1]
   local gamma = outFld[2]
   local gamma_inv = outFld[3]

   local localRange = p_over_gamma:localRange()
   ffiC.gkyl_calc_sr_vars_init_p_vars(self._velGrid._zero, self._velBasis._zero, 
      localRange, p_over_gamma._zero, gamma._zero, gamma_inv._zero)
end

-- Calculate Gamma^2 = 1/(1 - V^2/c^2) advance method.
function CalcSRVars:_advanceGamma2(tCurr, inFld, outFld)

   local V_drift = inFld[1]
   local Gamma2 = outFld[1]

   local localRange = Gamma2:localRange()
   ffiC.gkyl_calc_sr_vars_Gamma2(self._confBasis._zero, self._phaseBasis._zero, 
      localRange, V_drift._zero, Gamma2._zero)
end

-- Calculate Gamma = 1/sqrt(1 - V^2/c^2) advance method.
function CalcSRVars:_advanceGamma(tCurr, inFld, outFld)

   local V_drift = inFld[1]
   local Gamma = outFld[1]

   local localRange = Gamma:localRange()
   ffiC.gkyl_calc_sr_vars_Gamma(self._confBasis._zero, self._phaseBasis._zero, 
      localRange, V_drift._zero, Gamma._zero)
end

-- Calculate Gamma_inv = sqrt(1 - V^2/c^2) advance method.
function CalcSRVars:_advanceGammaInv(tCurr, inFld, outFld)

   local V_drift = inFld[1]
   local Gamma_inv = outFld[1]

   local localRange = Gamma_inv:localRange()
   ffiC.gkyl_calc_sr_vars_Gamma_inv(self._confBasis._zero, self._phaseBasis._zero, 
      localRange, V_drift._zero, Gamma_inv._zero)
end

function CalcSRVars:_advance(tCurr, inFld, outFld)
   -- init: Initialize SR variables (p/gamma, gamma, gamma_inv)
   -- Gamma2: Calculate Gamma^2 = 1/(1 - V^2/c^2)
   -- Gamma: Calculate Gamma = 1/sqrt(1 - V^2/c^2)
   -- Gamma_inv: Calculate 1/Gamma = sqrt(1 - V^2/c^2)
   self.advanceFunc(tCurr, inFld, outFld)
end

return CalcSRVars
