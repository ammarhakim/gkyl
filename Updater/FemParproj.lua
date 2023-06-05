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
local Grid        = require "Grid"
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
   self.xLCFS = tbl.xLCFS
   self.doubleBC = false

   if self._grid:ndim() ==3 and self._grid:numCells(2) == 1 and self.xLCFS then
      self.doubleBC = true
   end

   local weightFld  = tbl.weight
   local isWeighted = weightFld ~= nil
   if not isWeighted then  -- Create a dummy weight (not used).
      weightFld = DataStruct.Field {
         onGrid = self._grid,  numComponents = self._basis:numBasis(),
         ghost  = {1,1},       useDevice = false,
      }
      weightFld:clear(0.0)
   end


   if self.doubleBC then
      self._zero = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, false,
                                                 isWeighted, weightFld._zero, GKYL_USE_GPU or 0),
                       ffiC.gkyl_fem_parproj_release)

      self._zeroCore = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, true,
                                                 isWeighted, weightFld._zero, GKYL_USE_GPU or 0),
                       ffiC.gkyl_fem_parproj_release)
      self._advance = function(tCurr, inFld, outFld) return self:_avanceDoubleBC(tCurr,inFld, outFld) end
      self._advanceOnDevice = function(tCurr, inFld, outFld) return self:_advanceOnDeviceAxisymmetric(tCurr,inFld, outFld) end
      self.inSOL = DataStruct.Field { -- May need to add some arguments to this ield creation
         onGrid = self._grid,  numComponents = self._basis:numBasis(),
         ghost  = {1,1},       useDevice = false,
         metaData      = {polyOrder = self._basis:polyOrder(),
                          basisType = self._basis:id(),},
      }
      self.outSOL = DataStruct.Field { -- May need to add some arguments to this ield creation
         onGrid = self._grid,  numComponents = self._basis:numBasis(),
         ghost  = {1,1},       useDevice = false,
         metaData      = {polyOrder = self._basis:polyOrder(),
                          basisType = self._basis:id(),},
      }
      self:setSplitRanges(self.inSOL)
   else
      self._zero = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, self._isParPeriodic,
                                                 isWeighted, weightFld._zero, GKYL_USE_GPU or 0),
                       ffiC.gkyl_fem_parproj_release)
      self._advance = function(tCurr, inFld, outFld)
         return self:advanceSingleBC(tCurr,inFld, outFld)
      end
      self._advanceOnDevice = function(tCurr, inFld, outFld) return self:_advanceOnDeviceSingleBC(tCurr,inFld, outFld) end
   end

end

function FemParproj:setSplitRanges(sampleFld)
   local gridRange = self._grid:globalRange()
   local xLCFS=self.xLCFS
   local coordLCFS = {xLCFS-1.e-7}
   local idxLCFS    = {-9}
   local xGridIngr = self._grid:childGrid({1})
   local xGrid = Grid.RectCart {
      lower = xGridIngr.lower,  periodicDirs  = xGridIngr.periodicDirs,
      upper = xGridIngr.upper,  decomposition = xGridIngr.decomposition,
      cells = xGridIngr.cells,
   }
   xGrid:findCell(coordLCFS, idxLCFS)
   local globalCore = gridRange:shorten(1, idxLCFS[1])
   local globalSOLTemp = gridRange:shortenFromBelow(1, self._grid:numCells(1)-idxLCFS[1]+1)
   local lv, uv = globalSOLTemp:lowerAsVec(), globalSOLTemp:upperAsVec()
   self.globalSOL = gridRange:subRange(lv,uv)

   local localRange = sampleFld:localRange()
   local localSOLRangeTemp = sampleFld:localRange():intersect(self.globalSOL)
   local localCoreRangeTemp = sampleFld:localRange():intersect(globalCore)
   local lv, uv = localSOLRangeTemp:lowerAsVec(), localSOLRangeTemp:upperAsVec()
   self.localSOLRange = localRange:subRange(lv,uv)
   local lv, uv = localCoreRangeTemp:lowerAsVec(), localCoreRangeTemp:upperAsVec()
   self.localCoreRange = localRange:subRange(lv,uv)
end

function FemParproj:_avanceDoubleBC(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   self.inSOL:copy(rhoIn)
   self.inSOL:write('inSOL.bp', 0.0, 2, false)

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, self.inSOL._zero)
   ffiC.gkyl_fem_parproj_solve(self._zero, self.outSOL._zero)
   self.outSOL:write('outSOL.bp', 0.0, 2, false)

   rhoIn:periodicCopyInDir(3)
   rhoIn:write('rhoIn.bp', 0.0, 2, false)
   ffiC.gkyl_fem_parproj_set_rhs(self._zeroCore, rhoIn._zero)
   ffiC.gkyl_fem_parproj_solve(self._zeroCore, qOut._zero)
   qOut:write('qOut.bp', 0.0, 2, false)

   self.outSOL:write('outSOL2.bp', 0.0, 2, false)
   qOut:copyRangeToRange(self.outSOL,  self.localSOLRange, self.localSOLRange)
   qOut:write('qOut2.bp', 0.0, 2, false)
end


function FemParproj:advanceSingleBC(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   if self._isParPeriodic then rhoIn:sync(true) end

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, rhoIn._zero)
   ffiC.gkyl_fem_parproj_solve(self._zero, qOut._zero)
end


function FemParproj:_advanceOnDeviceSingleBC(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   if self._isParPeriodic then rhoIn:sync(true) end

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, rhoIn._zeroDevice)
   ffiC.gkyl_fem_parproj_solve(self._zero, qOut._zeroDevice)
end

function FemParproj:printDevDiagnostics() end

return FemParproj
