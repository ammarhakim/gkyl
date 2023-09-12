-- Gkyl ------------------------------------------------------------------------
--
-- Solve the perpendicular Poisson/Helmholtz equation
--   -nabla . (epsilon * nabla_perp phi) - kSq*phi = rho
-- using the continuous finite element method.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_fem_poisson_perp gkyl_fem_poisson_perp;

/**
 * Create new updater to solve the Helmholtz problem
 *   - nabla . (epsilon * nabla phi) - kSq * phi = rho
 * using a FEM to ensure phi is continuous. This solver is also
 * used as a Poisson solver by passing a zero kSq. The input is the
 * DG field rho, which is translated to FEM. The output is the
 * DG field phi, after we've translated the FEM solution to DG.
 * Free using gkyl_fem_poisson_perp_release method.
 *
 * @param solve_range Range in which to perform the projection operation.
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param bcs Boundary conditions.
 * @param epsilon Spatially varying permittivity tensor.
 * @param kSq Squared wave number (factor multiplying phi in Helmholtz eq).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_poisson_perp* gkyl_fem_poisson_perp_new(
  const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis basis, struct gkyl_poisson_bc *bcs, struct gkyl_array *epsilon,
  struct gkyl_array *kSq, bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM poisson updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_poisson_perp_set_rhs(gkyl_fem_poisson_perp* up, struct gkyl_array *rhsin);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_poisson_perp_solve(gkyl_fem_poisson_perp* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_perp_release(struct gkyl_fem_poisson_perp *up);
]]

-- Boundary condition updater.
local FemPoissonPerp = Proto(UpdaterBase)

function FemPoissonPerp:init(tbl)
   FemPoissonPerp.super.init(self, tbl) -- Setup base object.

   self._grid   = assert(tbl.onGrid, "Updater.FemPoissonPerp: Must specify grid to use with 'onGrid'.")
   self._basis  = assert(tbl.basis, "Updater.FemPoissonPerp: Must specify the basis in 'basis'.")
   local eps    = assert(tbl.epsilon, "Updater.FemPoissonPerp: Must specify the permittivity 'epsilon'.")
   local kSq    = tbl.kSq  -- Wave number squared.
   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local ndim = self._grid:ndim()

   local function translateBcType(bcTypeIn)
      -- These have to match gkyl_poisson_bc_type in gkylzero/zero/gkyl_fem_poisson_perp.h.
      local bcKey = {GKYL_POISSON_PERIODIC  = 0,
                     GKYL_POISSON_DIRICHLET = 1,
                     GKYL_POISSON_NEUMANN   = 2,
      }
          if bcTypeIn == "P" then return bcKey["GKYL_POISSON_PERIODIC"]
      elseif bcTypeIn == "D" then return bcKey["GKYL_POISSON_DIRICHLET"]
      elseif bcTypeIn == "N" then return bcKey["GKYL_POISSON_NEUMANN"]
      else assert(false, "Updater.FemPoissonPerp: boundary condition type must be specified by one of 'P', 'D', 'N'.")
      end
   end

   local bc_zero = ffi.new("struct gkyl_poisson_bc")
   if tbl.bcLower and tbl.bcUpper then
      for d = 1,ndim do
         bc_zero.lo_type[d-1] = translateBcType(tbl.bcLower[d].T)
         bc_zero.up_type[d-1] = translateBcType(tbl.bcUpper[d].T)

         if tbl.bcLower[d].T ~= "P" then
            bc_zero.lo_value[d-1].v[0] = tbl.bcLower[d].V
         end
         if tbl.bcUpper[d].T ~= "P" then
            bc_zero.up_value[d-1].v[0] = tbl.bcUpper[d].V
         end
      end
   else
      assert(false, "Updater.FemPoissonPerp: must specify 'bcLower' and 'bcUpper'.")
   end

   assert(type(eps) == 'table', "Updater.FemPoissonPerp: permittivity epsilon must be a CartField.")
   local eps_p = self._useGPU and eps._zeroDevice or eps._zero
   local kSq_p = nil
   if kSq then
      assert(type(kSq) == 'table', "Updater.FemPoissonPerp: squared wave-number kSq must be a CartField.")
      kSq_p = self._useGPU and kSq._zeroDevice or kSq._zero
   end

   local localRange = eps:localRange()

   self._zero = ffi.gc(ffiC.gkyl_fem_poisson_perp_new(localRange, self._grid._zero, self._basis._zero, bc_zero,
                                                      eps_p, kSq_p, self._useGPU),
                       ffiC.gkyl_fem_poisson_perp_release)
end

function FemPoissonPerp:_advance(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemPoissonPerp.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemPoissonPerp.advance: Must-specify an output field")

   ffiC.gkyl_fem_poisson_perp_set_rhs(self._zero, rhoIn._zero)
   ffiC.gkyl_fem_poisson_perp_solve(self._zero, qOut._zero)
end

function FemPoissonPerp:_advanceOnDevice(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemPoissonPerp.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemPoissonPerp.advance: Must-specify an output field")

   ffiC.gkyl_fem_poisson_perp_set_rhs(self._zero, rhoIn._zeroDevice)
   ffiC.gkyl_fem_poisson_perp_solve(self._zero, qOut._zeroDevice)
end

function FemPoissonPerp:printDevDiagnostics() end

return FemPoissonPerp
