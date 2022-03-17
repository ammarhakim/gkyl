--
-- Updater to compute RHS or forward Euler update for hyperbolic
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Basis = require "Basis.BasisCdef"
local DataStruct = require "DataStruct"
local EqBase = require "Eq.EqBase"
local Grid = require "Grid.RectCart"
local CartField = require "DataStruct.CartField"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local cuda = nil
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
   cuAlloc = require "Cuda.Alloc"
end

ffi.cdef [[ 
// Object type
struct gkyl_hyper_dg {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_basis; // number of basis functions
  int num_up_dirs; // number of update directions
  int update_dirs[7]; // directions to update
  // zero_flux_flags[d] == 1 means zero-flux BC in 'd'
  int zero_flux_flags[7];
  int update_vol_term; // should we update volume term?
  const struct gkyl_dg_eqn *equation; // equation object

  uint32_t flags;
  struct gkyl_hyper_dg *on_dev; // pointer to itself or device data
};
// Object type
typedef struct gkyl_hyper_dg gkyl_hyper_dg;

/**
 * Create new updater to update equations using DG algorithm.
 *
 * @param grid Grid object
 * @param basis Basis functions
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 * @param zero_flux_flags Flags to indicate if direction has zero-flux BCs
 * @param update_vol_term Set to 0 to skip volume update
 */
gkyl_hyper_dg* gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[7], int zero_flux_flags[7],
  int update_vol_term, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param hdg Hyper DG updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_advance(gkyl_hyper_dg *hdg, const struct gkyl_range *update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs);

void gkyl_hyper_dg_advance_cu(gkyl_hyper_dg* hdg, const struct gkyl_range *update_range,
  const struct gkyl_array* fIn, struct gkyl_array* cflrate,
  struct gkyl_array* rhs);

]]

-- Hyperbolic DG solver updater object
local HyperDisCont = Proto(UpdaterBase)

function HyperDisCont:init(tbl)
   HyperDisCont.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.HyperDisCont: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.HyperDisCont: Must specify basis functions to use using 'basis'")

   -- By default, clear output field before incrementing with vol/surf updates.
   self._clearOut = xsys.pickBool(tbl.clearOut, true)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._onGrid:ndim()

   -- Equation to solve.
   self._equation = assert(tbl.equation, "Updater.HyperDisCont: Must provide equation object using 'equation'")

   -- By default, update all directions.
   self._updateDirs = {}
   for d = 1, self._ndim do
      self._updateDirs[d] = d
   end
   -- Read in which directions we are to update.
   if tbl.updateDirections then
      self._updateDirs = tbl.updateDirections
   end

   -- Set zero flux direction flags.
   self._zeroFluxFlags = {}
   for d = 1, self._ndim do self._zeroFluxFlags[d] = false end
   if tbl.zeroFluxDirections then
      for i, d in ipairs(tbl.zeroFluxDirections) do self._zeroFluxFlags[d] = true end
   end

   -- Flag to turn on/off volume term.
   self._updateVolumeTerm = xsys.pickBool(tbl.updateVolumeTerm, true)

   -- Flag to indicate the use of local or global upwind fluxes. Local
   -- does not need to do a reduction of the maximum speed (saves an Mpi.Allreduce).
   self._globalUpwind = xsys.pickBool(tbl.globalUpwind, true)

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.HyperDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- Maximum characteristic velocities for use in penalty based fluxes.
   self._maxs, self._maxsOld, self._maxsLocal = Lin.Vec(self._ndim), Lin.Vec(self._ndim), Lin.Vec(self._ndim)
   for d = 1, self._ndim do
      -- Very first step the penalty term will be zero. However, in
      -- subsequent steps the maximum speed from the previous step
      -- will be used.
      self._maxs[d] = 0.0
   end
   self._noPenaltyFlux = xsys.pickBool(tbl.noPenaltyFlux, false)

   self._isFirst = true
   self._auxFields = {} -- Auxilliary fields passed to eqn object.
   self._perpRangeDecomp = {} -- Perp ranges in each direction.
   
   local upd = ffi.new("int[6]", self._updateDirs)
   for i = 1, #self._updateDirs do
      upd[i-1] = self._updateDirs[i]-1
   end
   local zf = ffi.new("int[6]", self._zeroFluxFlags)
   if self._equation._zero then 
     self._zero = ffiC.gkyl_hyper_dg_new(self._onGrid._zero, self._basis._zero, self._equation._zero, #self._updateDirs, upd, zf, self._updateVolumeTerm, GKYL_USE_GPU or 0)
   end

   return self
end

-- advance method
function HyperDisCont:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local qIn = assert(inFld[1], "HyperDisCont.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "HyperDisCont.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "HyperDisCont.advance: Must pass cflRate field in output table")

   -- pass aux fields to equation object
   for i = 1, #inFld-1 do
      self._auxFields[i] = inFld[i+1]
   end
   self._equation:setAuxFields(self._auxFields)

   local ndim = grid:ndim()

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local localRange = qRhsOut:localRange()
   local globalRange = qRhsOut:globalRange()
   local qInIdxr, qRhsOutIdxr = qIn:genIndexer(), qRhsOut:genIndexer() -- indexer functions into fields
   local cflRateByCellIdxr = cflRateByCell:genIndexer()

   if self._zero then
      if self._clearOut then qRhsOut:clear(0.0) end
      ffiC.gkyl_hyper_dg_advance(self._zero, localRange, qIn._zero.arr, cflRateByCell._zero.arr, qRhsOut._zero.arr)
   else
      -- to store grid info
      local dxp, dxm = Lin.Vec(ndim), Lin.Vec(ndim) -- cell shape on right/left
      local xcp, xcm = Lin.Vec(ndim), Lin.Vec(ndim) -- cell center on right/left
      local idxp, idxm = Lin.IntVec(ndim), Lin.IntVec(ndim) -- index on right/left

      -- pointers for (re)use in update
      local qInM, qInP = qIn:get(1), qIn:get(1)
      local qRhsOutM, qRhsOutP = qRhsOut:get(1), qRhsOut:get(1)
      local cflRateByCellP = cflRateByCell:get(1)
      local cflRateByCellM = cflRateByCell:get(1)

      -- This flag is needed as the volume integral already contains
      -- contributions from all directions. Hence, we must only
      -- accumulate the volume contribution once, skipping it for other
      -- directions.
      local firstDir = true

      -- Use maximum characteristic speeds from previous step as penalty.
      for d = 1, ndim do
         if self._noPenaltyFlux then 
            self._maxsOld[d] = 0.0
         else
            self._maxsOld[d] = self._maxs[d]
         end
         self._maxsLocal[d] = 0.0 -- Reset to get new values in this step.
      end

      local tId = grid:subGridSharedId() -- Local thread ID.

      -- Clear output field before computing vol/surf increments.
      if self._clearOut then qRhsOut:clear(0.0) end
      -- Accumulate contributions from volume and surface integrals.
      local cflRate
      -- Iterate through updateDirs backwards so that a zero flux dir is first in kinetics.
      for i = #self._updateDirs, 1, -1 do 
         local dir = self._updateDirs[i]
         -- Lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell).
         local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
         local dirLoSurfIdx, dirUpSurfIdx = dirLoIdx, dirUpIdx
         
         -- Compute loop bounds for zero flux direction.
         if self._zeroFluxFlags[dir] then
            local dirGlobalLoIdx, dirGlobalUpIdx = globalRange:lower(dir), globalRange:upper(dir)+1
            if dirLoIdx == dirGlobalLoIdx then 
               dirLoSurfIdx = dirLoIdx+1
            end
            if dirUpIdx == dirGlobalUpIdx then 
               dirUpSurfIdx = dirUpIdx-1
            end
         end

         if self._isFirst then
            self._perpRangeDecomp[dir] = LinearDecomp.LinearDecompRange {
               range = localRange:shorten(dir), -- range orthogonal to 'dir'
               numSplit = grid:numSharedProcs(),
               threadComm = self:getSharedComm()
            }
         end
         local perpRangeDecomp = self._perpRangeDecomp[dir]

         -- outer loop is over directions orthogonal to 'dir' and inner
         -- loop is over 1D slice in `dir`.
         for idx in perpRangeDecomp:rowMajorIter(tId) do
            idx:copyInto(idxp); idx:copyInto(idxm)

            for i = dirLoIdx, dirUpIdx do -- this loop is over edges
               idxm[dir], idxp[dir] = i-1, i -- cell left/right of edge 'i'

               grid:setIndex(idxm)
               grid:getDx(dxm)
               grid:cellCenter(xcm)

               grid:setIndex(idxp)
               grid:getDx(dxp)
               grid:cellCenter(xcp)

               qIn:fill(qInIdxr(idxm), qInM)
               qIn:fill(qInIdxr(idxp), qInP)

               qRhsOut:fill(qRhsOutIdxr(idxm), qRhsOutM)
               qRhsOut:fill(qRhsOutIdxr(idxp), qRhsOutP)
               cflRateByCell:fill(cflRateByCellIdxr(idxm), cflRateByCellM)
               cflRateByCell:fill(cflRateByCellIdxr(idxp), cflRateByCellP)

               if firstDir and i<=dirUpIdx-1 and self._updateVolumeTerm then
                  cflRate = self._equation:volTerm(xcp, dxp, idxp, qInP, qRhsOutP)
                  cflRateByCellP:data()[0] = cflRateByCellP:data()[0] + cflRate
               end
               if i >= dirLoSurfIdx and i <= dirUpSurfIdx then
                  local maxs = self._equation:surfTerm(
           	  dir, 1.0, 1.0, xcm, xcp, dxm, dxp, self._maxsOld[dir], idxm, idxp, qInM, qInP, qRhsOutM, qRhsOutP)
                  self._maxsLocal[dir] = math.max(self._maxsLocal[dir], maxs)
               else
                  if self._zeroFluxFlags[dir] then
                     -- we need to give equations a chance to apply partial
                     -- surface updates even when the zeroFlux BCs have been
                     -- applied
                     self._equation:boundarySurfTerm(
           	     dir, xcm, xcp, dxm, dxp, self._maxsOld[dir], idxm, idxp, qInM, qInP, qRhsOutM, qRhsOutP)
                  end
               end
            end
         end
         if firstDir then cflRateByCell:sync() end
         firstDir = false
      end

      -- Determine largest amax across processors.
      if self._globalUpwind then
         Mpi.Allreduce(
            self._maxsLocal:data(), self._maxs:data(), ndim, Mpi.DOUBLE, Mpi.MAX, self:getComm())
      end
      self._isFirst = false
   end
end

function HyperDisCont:_advanceOnDevice(tCurr, inFld, outFld)
   local grid = self._onGrid

   local qIn = assert(inFld[1], "HyperDisCont.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "HyperDisCont.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "HyperDisCont.advance: Must pass cflRate field in output table")

   -- pass aux fields to equation object
   for i = 1, #inFld-1 do
      self._auxFields[i] = inFld[i+1]
   end
   self._equation:setAuxFieldsOnDevice(self._auxFields)

   local localRange = qRhsOut:localRange()

   qRhsOut:clear(0.0)
   ffiC.gkyl_hyper_dg_advance_cu(self._zero, localRange, qIn._zeroDevice.arr, cflRateByCell._zeroDevice.arr, qRhsOut._zeroDevice.arr)
end

-- set up pointers to dt and cflRateByCell
function HyperDisCont:setDtAndCflRate(dt, cflRateByCell)
   HyperDisCont.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

return HyperDisCont
