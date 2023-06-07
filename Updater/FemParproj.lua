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

   local function createField(grid, basis, ghostCells, useDevice)
      local useDevice = xsys.pickBool(useDevice, GKYL_USE_GPU)
      local fld = DataStruct.Field {
         onGrid           = grid,
         numComponents    = basis:numBasis(),
         ghost            = ghostCells,
         metaData         = {polyOrder = basis:polyOrder(),
                             basisType = basis:id()},
         useDevice = useDevice,
      }
      fld:clear(0.0)
      return fld
   end

   local weightFld  = tbl.weight
   local isWeighted = weightFld ~= nil
   if not isWeighted then  -- Create a dummy weight (not used).
      weightFld = createField(self._grid, self._basis, {1,1}, false)
   end


   if self.doubleBC then
      self._zero = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, false,
                                                 isWeighted, weightFld._zero, GKYL_USE_GPU or 0),
                       ffiC.gkyl_fem_parproj_release)

      self._zeroCore = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, true,
                                                 isWeighted, weightFld._zero, GKYL_USE_GPU or 0),
                       ffiC.gkyl_fem_parproj_release)
      self.advanceFunc = function(tCurr, inFld, outFld) self:_advanceDoubleBC(tCurr,inFld, outFld) end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld) self:_advanceOnDeviceDoubleBC(tCurr,inFld, outFld) end
      self.inFldRegion2 = createField(self._grid, self._basis, {1,1})
      self.outFldRegion2 = createField(self._grid, self._basis, {1,1})
      self:setSplitRanges(self.inFldRegion2)
   else
      self._zero = ffi.gc(ffiC.gkyl_fem_parproj_new(self._grid._zero, self._basis._zero, self._isParPeriodic,
                                                 isWeighted, weightFld._zero, GKYL_USE_GPU or 0),
                       ffiC.gkyl_fem_parproj_release)
      self.advanceFunc = function(tCurr, inFld, outFld) self:advanceSingleBC(tCurr,inFld, outFld) end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld) self:_advanceOnDeviceSingleBC(tCurr,inFld, outFld) end
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
   local globalSOL = gridRange:subRange(lv,uv)

   local localRange = sampleFld:localRange()
   local localSOLRangeTemp = sampleFld:localRange():intersect(globalSOL)
   local localCoreRangeTemp = sampleFld:localRange():intersect(globalCore)
   local lv, uv = localSOLRangeTemp:lowerAsVec(), localSOLRangeTemp:upperAsVec()
   self.localRegion2Range = localRange:subRange(lv,uv)
   local lv, uv = localCoreRangeTemp:lowerAsVec(), localCoreRangeTemp:upperAsVec()
   self.localRegion1Range = localRange:subRange(lv,uv)
end

function FemParproj:_advanceDoubleBC(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   self.inFldRegion2:copy(rhoIn)
   local fldGrid = rhoIn:grid()
   local periodicDirs = fldGrid:getPeriodicDirs()
   local modifiedPeriodicDirs = {3}
   fldGrid:setPeriodicDirs(modifiedPeriodicDirs)
   rhoIn:sync()

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, self.inFldRegion2._zero)
   ffiC.gkyl_fem_parproj_solve(self._zero, self.outFldRegion2._zero)

   ffiC.gkyl_fem_parproj_set_rhs(self._zeroCore, rhoIn._zero)
   ffiC.gkyl_fem_parproj_solve(self._zeroCore, qOut._zero)

   qOut:copyRange(self.outFldRegion2, self.localRegion2Range)
   fldGrid:setPeriodicDirs(periodicDirs)
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


function FemParproj:_advanceOnDeviceDoubleBC(tCurr, inFld, outFld)
   local rhoIn = assert(inFld[1], "FemParproj.advance: Must-specify an input field")
   local qOut  = assert(outFld[1], "FemParproj.advance: Must-specify an output field")

   self.inFldRegion2:copy(rhoIn)
   local fldGrid = rhoIn:grid()
   local periodicDirs = fldGrid:getPeriodicDirs()
   local modifiedPeriodicDirs = {3}
   fldGrid:setPeriodicDirs(modifiedPeriodicDirs)
   rhoIn:sync()

   ffiC.gkyl_fem_parproj_set_rhs(self._zero, self.inFldRegion2._zeroDevice)
   ffiC.gkyl_fem_parproj_solve(self._zero, self.outFldRegion2._zeroDevice)

   ffiC.gkyl_fem_parproj_set_rhs(self._zeroCore, rhoIn._zeroDevice)
   ffiC.gkyl_fem_parproj_solve(self._zeroCore, qOut._zeroDevice)

   qOut:copyRange(self.outFldRegion2, self.localRegion2Range)
   fldGrid:setPeriodicDirs(periodicDirs)
end

function FemParproj:_advance(tCurr, inFld, outFld)
   self.advanceFunc(tCurr, inFld, outFld)
end

function FemParproj:_advanceOnDevice(tCurr, inFld, outFld)
   self.advanceOnDeviceFunc(tCurr, inFld, outFld)
end

function FemParproj:printDevDiagnostics() end

return FemParproj
