-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute RHS or forward Euler update for hyperbolic
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local xsys = require "xsys"

local UPDATE_BOTH  = 0
local UPDATE_LEFT  = 1
local UPDATE_RIGHT = 2

-- Hyperbolic DG solver updater object
local HyperDisContCellBased = Proto(UpdaterBase)

function HyperDisContCellBased:init(tbl)
   HyperDisContCellBased.super.init(self, tbl)

   -- read data from input file
   self._onGrid = assert(tbl.onGrid, "Updater.HyperDisContCellBased: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.HyperDisContCellBased: Must specify basis functions to use using 'basis'")

   -- by default, clear output field before incrementing with vol/surf updates
   self._clearOut = xsys.pickBool(tbl.clearOut, true)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._onGrid:ndim()

   -- equation to solve
   self._equation = assert(tbl.equation, "Updater.HyperDisContCellBased: Must provide equation object using 'equation'")

   -- by default, update all directions
   self._updateDirs = {}
   for d = 1, self._ndim do
      self._updateDirs[d] = d
   end
   -- read in which directions we are to update
   if tbl.updateDirections then
      self._updateDirs = tbl.updateDirections
   end

   -- set zero flux direction flags
   self._zeroFluxFlags = {}
   for d = 1, self._ndim do
      self._zeroFluxFlags[d] = false
   end
   if tbl.zeroFluxDirections then
      for i, d in ipairs(tbl.zeroFluxDirections) do
         self._zeroFluxFlags[d] = true
      end
   end

   -- flag to turn on/off volume term
   self._updateVolumeTerm = xsys.pickBool(tbl.updateVolumeTerm, true)

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.HyperDisContCellBased: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- maximum characteristic velocities for use in penalty based
   -- fluxes
   self._maxs, self._maxsOld, self._maxsLocal = Lin.Vec(self._ndim), Lin.Vec(self._ndim), Lin.Vec(self._ndim)
   for d = 1, self._ndim do
      -- very first step the penalty term will be zero. However, in
      -- subsequent steps the maximum speed from the previous step
      -- will be used
      self._maxs[d] = 0.0
   end

   self._isFirst = true
   self._auxFields = {} -- auxilliary fields passed to eqn object
   self._perpRangeDecomp = {} -- perp ranges in each direction      

   return self
end

-- advance method
function HyperDisContCellBased:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local cflRateByCell = self._cflRateByCell

   local qIn = assert(inFld[1], "HyperDisContCellBased.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "HyperDisContCellBased.advance: Must specify an output field")

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

   -- to store grid info
   local dxC, dxL, dxR = Lin.Vec(ndim), Lin.Vec(ndim), Lin.Vec(ndim) -- cell shape on right/left
   local xcC, xcL, xcR = Lin.Vec(ndim), Lin.Vec(ndim), Lin.Vec(ndim) -- cell center on right/left
   local idxC, idxL, idxR = Lin.IntVec(ndim), Lin.IntVec(ndim), Lin.IntVec(ndim) -- index on right/left

   -- pointers for (re)use in update
   local qInC = qIn:get(1)
   local qInL, qInR = qIn:get(1), qIn:get(1)
   local qRhsOutC = qRhsOut:get(1)
   local qRhsOutL, qRhsOutR = qRhsOut:get(1), qRhsOut:get(1)
   local cflRateByCellC = cflRateByCell:get(1)
   local cflRateByCellR = cflRateByCell:get(1)
   local cflRateByCellL = cflRateByCell:get(1)

   -- This flag is needed as the volume integral already contains
   -- contributions from all directions. Hence, we must only
   -- accumulate the volume contribution once, skipping it for other
   -- directions
   local firstDir = true

   -- use maximum characteristic speeds from previous step as penalty
   for d = 1, ndim do
      self._maxsOld[d] = self._maxs[d]
      self._maxsLocal[d] = 0.0 -- reset to get new values in this step
   end

   local tId = grid:subGridSharedId() -- local thread ID

   -- clear output field before computing vol/surf increments
   if self._clearOut then qRhsOut:clear(0.0) end
   -- accumulate contributions from volume and surface integrals
   local cflRate

   local rangeDecomp = LinearDecomp.LinearDecompRange {
      range = globalRange, numSplit = grid:numSharedProcs() }


   -- cell-based loop
   for idxC in rangeDecomp:rowMajorIter(tId) do

      -- volume update
      grid:setIndex(idxC)
      grid:getDx(dxC)
      grid:cellCenter(xcC)

      qIn:fill(qInIdxr(idxC), qInC)
      
      qRhsOut:fill(qRhsOutIdxr(idxC), qRhsOutC)
      cflRateByCell:fill(cflRateByCellIdxr(idxC), cflRateByCellC)

      cflRate = self._equation:volTerm(xcC, dxC, idxC, qInC, qRhsOutC)
      cflRateByCellC:data()[0] = cflRateByCellC:data()[0] + cflRate

      -- surface update
      for i = 1, self._updateDirs do
         local dir = self._updateDirs[i]
          
         -- get indexes of cells to left (L) and right (R)
         idxL[dir] = idxC[dir] - 1
         idxR[dir] = idxC[dir] + 1

         grid:setIndex(idxL)
         grid:getDx(dxL)
         grid:cellCenter(xcL)
         
         grid:setIndex(idxR)
         grid:getDx(dxR)
         grid:cellCenter(xcR)
         
         -- get L and R ptrs
         qIn:fill(qInIdxr(idxL), qInL)
         qIn:fill(qInIdxr(idxR), qInR)
         qRhsOut:fill(qRhsOutIdxr(idxL), qRhsOutL)
         qRhsOut:fill(qRhsOutIdxr(idxR), qRhsOutR)
         cflRateByCell:fill(cflRateByCellIdxr(idxL), cflRateByCellL)
         cflRateByCell:fill(cflRateByCellIdxr(idxR), cflRateByCellR)

         local cflL = cflRateByCellL:data()[0]*dt/.9 -- .9 here is conservative, but we are using dt from the prev step
         local cflR = cflRateByCellR:data()[0]*dt/.9 -- .9 here is conservative, but we are using dt from the prev step

         -- left surface update
         if not (self._zeroFluxFlags[dir] and idxC[dir] == globalRange:lower(dir)) then
	    local maxs = self._equation:surfTerm(
	       dir, cflL, cflC, xcL, xcC, dxL, dxC, self._maxsOld[dir], idxL, idxC, qInL, qInC, qRhsOutL, qRhsOutC, UPDATE_RIGHT)
	    self._maxsLocal[dir] = math.max(self._maxsLocal[dir], maxs)
         else
            if self._zeroFluxFlags[dir] then
               -- we need to give equations a chance to apply partial
               -- surface updates even when the zeroFlux BCs have been
               -- applied
               self._equation:boundarySurfTerm(
                  dir, xcL, xcC, dxL, dxC, self._maxsOld[dir], idxL, idxC, qInL, qInC, qRhsOutL, qRhsOutC, UPDATE_RIGHT)
            end
         end

         -- right surface update
         if not (self._zeroFluxFlags[dir] and idxC[dir] == globalRange:upper(dir)) then
	    local maxs = self._equation:surfTerm(
	       dir, cflC, cflR, xcC, xcR, dxC, dxR, self._maxsOld[dir], idxC, idxR, qInC, qInR, qRhsOutC, qRhsOutR, UPDATE_LEFT)
	    self._maxsLocal[dir] = math.max(self._maxsLocal[dir], maxs)
         else
            if self._zeroFluxFlags[dir] then
               -- we need to give equations a chance to apply partial
               -- surface updates even when the zeroFlux BCs have been
               -- applied
               self._equation:boundarySurfTerm(
                  dir, xcC, xcR, dxC, dxR, self._maxsOld[dir], idxC, idxR, qInC, qInR, qRhsOutC, qRhsOutR, UPDATE_LEFT)
            end
         end
      end
   end

   -- determine largest amax across processors
   Mpi.Allreduce(
      self._maxsLocal:data(), self._maxs:data(), ndim, Mpi.DOUBLE, Mpi.MAX, self:getComm())

   self._isFirst = false
end

return HyperDisContCellBased
