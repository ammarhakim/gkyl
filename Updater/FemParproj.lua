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
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param isparperiodic boolean indicating if parallel direction is periodic.
 * @param isweighted boolean indicating if wgt\=1.
 * @param weight multiplicative weight on left-side of the operator.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_fem_parproj* gkyl_fem_parproj_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis basis,
  bool isparperiodic, bool isweighted, const struct gkyl_array *weight,
  bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_parproj_set_rhs(gkyl_fem_parproj* up, const struct gkyl_array *rhsin);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve(gkyl_fem_parproj* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_release(gkyl_fem_parproj *up);
]]

-- Boundary condition updater.
local FemParproj = Proto(UpdaterBase)

function FemParproj:init(tbl)
   FemParproj.super.init(self, tbl) -- Setup base object.

   self._grid  = assert(tbl.onGrid, "Updater.FemParproj: Must specify grid to use with 'onGrid'.")
   self._basis = assert(tbl.basis, "Updater.FemParproj: Must specify the basis in 'basis'.")

   self._isParPeriodic = tbl.periodicParallelDir
   assert(self._isParPeriodic ~= nil, "Updater.FemParproj: Must specify if parallel direction is periodic with 'periodicParallelDir'.")
   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local weightFld  = tbl.weight
   local isWeighted = weightFld ~= nil
   if not isWeighted then  -- Create a dummy weight (not used).
      weightFld = DataStruct.Field {
         onGrid = self._grid,  numComponents = self._basis:numBasis(),
         ghost  = {1,1},       useDevice = false,
      }
      weightFld:clear(0.0)
   end

   self._zero = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, self._isParPeriodic,
                                                 isWeighted, weightFld._zero, self._useGPU),
                       ffiC.gkyl_fem_parproj_release)
end

function FemParproj:_advance(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   if self._isParPeriodic then rhoIn:sync(true) end

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, rhoIn._zero)
   ffiC.gkyl_fem_parproj_solve(self._zero, qOut._zero)
end

function FemParproj:_advanceOnDevice(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   if self._isParPeriodic then rhoIn:sync(true) end

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, rhoIn._zeroDevice)
   ffiC.gkyl_fem_parproj_solve(self._zero, qOut._zeroDevice)
end

function FemParproj:printDevDiagnostics() end

return FemParproj
