-- Gkyl ------------------------------------------------------------------------
--
-- Apply the FEM parallel projection operator to a DG field.
--
-- Primarily intended as a field solver in cdim=1 gyrokinetics, or as a
-- parallel smoother in cdim=3 gyrokinetics.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local DataStruct  = require "DataStruct"
local ffi         = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_fem_parproj gkyl_fem_parproj;

// Boundary condition types.
enum gkyl_fem_parproj_bc_type {
  GKYL_FEM_PARPROJ_PERIODIC = 0,
  GKYL_FEM_PARPROJ_DIRICHLET, // sets the value.
  GKYL_FEM_PARPROJ_NONE,      // does not enforce a BC.
};

/**
 * Create new updater to project a DG field onto the FEM (nodal) basis
 * in order to make the field continuous or, thanks to the option to pass
 * a multiplicative weight, solve 1D algebraic equations in which the output
 * field is continuous (but the input may not be). That is, we solve
 *    wgt*phi_{fem} \doteq rho_{dg}
 * where wgt is the weight field, phi_{fem} is the (continuous field)
 * we wish to compute, rho_{dg} is the (discontinuous) input source field,
 * and \doteq implies weak equality with respect to the FEM basis.
 * Free using gkyl_fem_parproj_release method.
 *
 * @param solve_range Range in which to perform the projection operation.
 * @param solve_range_ext solve_range with ghost cells (in z primarily) used
 *                        for the Dirichlet BC case.
 * @param basis Basis functions of the DG field.
 * @param bctype Type of boundary condition (see gkyl_fem_parproj_bc_type).
 * @param weight multiplicative weight on left-side of the operator.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_parproj* gkyl_fem_parproj_new(
  const struct gkyl_range *solve_range, const struct gkyl_range *solve_range_ext,
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight, bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_set_rhs(struct gkyl_fem_parproj* up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve(struct gkyl_fem_parproj* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_release(struct gkyl_fem_parproj *up);
]]

-- Boundary condition updater.
local FemParproj = Proto(UpdaterBase)

function FemParproj:init(tbl)
   FemParproj.super.init(self, tbl) -- Setup base object.

   self._grid  = assert(tbl.onGrid, "Updater.FemParproj: Must specify grid to use with 'onGrid'.")
   self._basis = assert(tbl.basis, "Updater.FemParproj: Must specify the basis in 'basis'.")
   local onField = assert(tbl.onField, "Updater.FemParproj: Must specify a sample field (to get ranges) in 'onField'.")

   local isParPeriodic = tbl.periodicParallelDir
   assert(isParPeriodic ~= nil, "Updater.FemParproj: Must specify if parallel direction is periodic with 'periodicParallelDir'.")
   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local weightFld = tbl.weight and tbl.weight._zero or nil

   assert((not (tbl.onRange and tbl.onExtRange)) or (tbl.onRange and tbl.onExtRange), "Updater.FemParproj: must specify both 'onRange' and 'onExtRange' or neither.")
   local solveRange    = tbl.onRange or onField:localRange()
   local solveExtRange = tbl.onExtRange or onField:localExtRange()

   local parBC = isParPeriodic and 0 or 2  -- MF 2023/08/01: Dirichlet not yet implemented.

   self._zero = ffi.gc(ffiC.gkyl_fem_parproj_new(solveRange, solveExtRange, self._basis._zero, parBC, weightFld, self._useGPU),
                       ffiC.gkyl_fem_parproj_release)
end

function FemParproj:_advance(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, rhoIn._zero, nil)
   ffiC.gkyl_fem_parproj_solve(self._zero, qOut._zero)
end

function FemParproj:_advanceOnDevice(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, rhoIn._zeroDevice, nil)
   ffiC.gkyl_fem_parproj_solve(self._zero, qOut._zeroDevice)
end

function FemParproj:printDevDiagnostics() end

return FemParproj
