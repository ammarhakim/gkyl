-- Gkyl ------------------------------------------------------------------------
--
-- Solve the Poisson equation -epsilon_0*del^2(phi)=rho using the
-- continuous finite element method.
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
typedef struct gkyl_fem_poisson gkyl_fem_poisson;

// Boundary condition types.
enum gkyl_poisson_bc_type {
  GKYL_POISSON_PERIODIC=0,
  GKYL_POISSON_DIRICHLET,  // sets the value.
  GKYL_POISSON_NEUMANN,  // sets the slope normal to the boundary.
  GKYL_POISSON_ROBIN,  // a combination of dirichlet and neumann.
};

// Boundary condition values. Dirichlet and Neumann use only one value,
// Robin uses 3, and periodic ignores the value.
struct gkyl_poisson_bc_value { double v[3]; };

struct gkyl_poisson_bc {
  enum gkyl_poisson_bc_type lo_type[3], up_type[3];
  struct gkyl_poisson_bc_value lo_value[3], up_value[3];
};

/**
 * Create new updater to solve the Helmholtz problem
 *   - nabla . (epsilon * nabla phi) - kSq * phi = rho
 * using a FEM to ensure phi is continuous. This solver is also
 * used as a Poisson solver by passing a zero kSq. The input is the
 * DG field rho, which is translated to FEM. The output is the
 * DG field phi, after we've translated the FEM solution to DG.
 * Free using gkyl_fem_poisson_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param bcs Boundary conditions.
 * @param epsilon_const Constant scalar value of the permittivity.
 * @param epsilon_var Spatially varying permittivity tensor.
 * @param kSq Squared wave number (factor multiplying phi in Helmholtz eq).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_fem_poisson* gkyl_fem_poisson_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis basis, struct gkyl_poisson_bc *bcs,
  double epsilon_const, struct gkyl_array *epsilon_var, struct gkyl_array *kSq, bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM poisson updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_poisson_set_rhs(gkyl_fem_poisson* up, struct gkyl_array *rhsin);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_poisson_solve(gkyl_fem_poisson* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_release(gkyl_fem_poisson *up);
]]

-- Boundary condition updater.
local FemPoisson = Proto(UpdaterBase)

function FemPoisson:init(tbl)
   FemPoisson.super.init(self, tbl) -- Setup base object.

   self._grid   = assert(tbl.onGrid, "Updater.FemPoisson: Must specify grid to use with 'onGrid'.")
   self._basis  = assert(tbl.basis, "Updater.FemPoisson: Must specify the basis in 'basis'.")
   local eps0   = assert(tbl.epsilon_0, "Updater.FemPoisson: Must specify the permittivity of space 'epsilon_0'.")
   local kSq    = tbl.kSq -- Wave number squared.
   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local ndim = self._grid:ndim()

   local function translateBcType(bcTypeIn)
      -- These have to match gkyl_poisson_bc_type in gkylzero/zero/gkyl_fem_poisson.h.
      local bcKey = {GKYL_POISSON_PERIODIC  = 0,
                     GKYL_POISSON_DIRICHLET = 1,
                     GKYL_POISSON_NEUMANN   = 2,
                     GKYL_POISSON_ROBIN     = 3,
      }
          if bcTypeIn == "P" then return bcKey["GKYL_POISSON_PERIODIC"]
      elseif bcTypeIn == "D" then return bcKey["GKYL_POISSON_DIRICHLET"]
      elseif bcTypeIn == "N" then return bcKey["GKYL_POISSON_NEUMANN"]
      elseif bcTypeIn == "R" then return bcKey["GKYL_POISSON_ROBIN"]
      else assert(false, "Updater.FemPoisson: boundary condition type must be specified by one of 'P', 'D', 'N', or 'R'.")
      end
   end

   local bc_zero = ffi.new("struct gkyl_poisson_bc")
   if tbl.bcLower and tbl.bcUpper then
      for d = 1,ndim do
         bc_zero.lo_type[d-1] = translateBcType(tbl.bcLower[d].T)
         bc_zero.up_type[d-1] = translateBcType(tbl.bcUpper[d].T)
         -- Robin BCs require 3 boundary values.
         if tbl.bcLower[d].T == "R" then
            for i = 1,3 do bc_zero.lo_value[d-1].v[i-1] = tbl.bcLower[d].V[i] end
         elseif tbl.bcLower[d].T ~= "P" then
            bc_zero.lo_value[d-1].v[0] = tbl.bcLower[d].V
         end
         if tbl.bcUpper[d].T == "R" then
            for i = 1,3 do bc_zero.up_value[d-1].v[i-1] = tbl.bcUpper[d].V[i] end
         elseif tbl.bcUpper[d].T ~= "P" then
            bc_zero.up_value[d-1].v[0] = tbl.bcUpper[d].V
         end
      end
   else
      assert(false, "Updater.FemPoisson: must specify 'bcLower' and 'bcUpper'.")
   end

   local eps0_const, eps0_var = nil, nil
   if type(eps0) == 'number' then
      eps0_const = eps0
   elseif type(eps0) == 'table' then
      -- eps0 must be an array with the tensor permittivity.
      eps0_var = self._useGPU and eps0._zeroDevice or eps0._zero
   end

   local kSq_p = nil
   if kSq then
      assert(type(kSq) == 'table', "Updater.FemPoissonPerp: squared wave-number kSq must be a CartField.")
      kSq_p = kSq._zero
   end

   self._zero = ffi.gc(ffiC.gkyl_fem_poisson_new(self._grid._zero, self._basis._zero, bc_zero,
                                                 eps0_const, eps0_var, kSq_p, self._useGPU),
                       ffiC.gkyl_fem_poisson_release)
end

function FemPoisson:_advance(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemPoisson.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemPoisson.advance: Must-specify an output field")

   ffiC.gkyl_fem_poisson_set_rhs(self._zero, rhoIn._zero)
   ffiC.gkyl_fem_poisson_solve(self._zero, qOut._zero)
end

function FemPoisson:_advanceOnDevice(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemPoisson.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemPoisson.advance: Must-specify an output field")

   ffiC.gkyl_fem_poisson_set_rhs(self._zero, rhoIn._zeroDevice)
   ffiC.gkyl_fem_poisson_solve(self._zero, qOut._zeroDevice)
end

function FemPoisson:printDevDiagnostics() end

return FemPoisson
