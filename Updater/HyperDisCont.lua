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

-- Hyperbolic DG solver updater object
local HyperDisCont = Proto(UpdaterBase)

function HyperDisCont:init(tbl)
   HyperDisCont.super.init(self, tbl)

   -- read data from input file
   self._onGrid = assert(tbl.onGrid, "Updater.HyperDisCont: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.HyperDisCont: Must specify basis functions to use using 'basis'")

   -- by default, compute forward Euler: if onlyIncrement is true,
   -- then only increments are computed. NOTE: The increments are NOT
   -- multiplied by dt.
   self._onlyIncrement = xsys.pickBool(tbl.onlyIncrement, false)

   -- by default, clear output field before incrementing with vol/surf updates
   self._clearOut = xsys.pickBool(tbl.clearOut, true)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._onGrid:ndim()

   -- equation to solve
   self._equation = assert(tbl.equation, "Updater.HyperDisCont: Must provide equation object using 'equation'")

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

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.HyperDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- maximum characteristic velocities for use in pentalty based
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
function HyperDisCont:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid

   local qIn = assert(inFld[1], "HyperDisCont.advance: Must specify an input field")
   local qOut = assert(outFld[1], "HyperDisCont.advance: Must specify an output field")

   -- pass aux fields to equation object
   for i = 1, #inFld-1 do
      self._auxFields[i] = inFld[i+1]
   end
   self._equation:setAuxFields(self._auxFields)

   local ndim = grid:ndim()

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local localRange = qOut:localRange()
   local globalRange = qOut:globalRange()
   local qInIdxr, qOutIdxr = qIn:genIndexer(), qOut:genIndexer() -- indexer functions into fields

   -- to store grid info
   local dxp, dxm = Lin.Vec(ndim), Lin.Vec(ndim) -- cell shape on right/left
   local xcp, xcm = Lin.Vec(ndim), Lin.Vec(ndim) -- cell center on right/left
   local idxp, idxm = Lin.IntVec(ndim), Lin.IntVec(ndim) -- index on right/left

   -- pointers for (re)use in update
   local qInM, qInP = qIn:get(1), qIn:get(1)
   local qOutM, qOutP = qOut:get(1), qOut:get(1)

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

   if self._clearOut then
      qOut:clear(0.0) -- clear output field before computing vol/surf increments
   end
   -- accumulate contributions from volume and surface integrals
  for _, dir in ipairs(self._updateDirs) do
      -- lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell)
      local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
      local dirLoSurfIdx, dirUpSurfIdx = dirLoIdx, dirUpIdx
      
      -- compute loop bounds for zero flux direction 
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
	    numSplit = grid:numSharedProcs() }
      end
      local perpRangeDecomp = self._perpRangeDecomp[dir]
      local tId = grid:subGridSharedId() -- local thread ID

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`.
      for idx in perpRangeDecomp:colMajorIter(tId) do
	 idx:copyInto(idxp); idx:copyInto(idxm)

   	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'

	    grid:setIndex(idxm)
	    for d = 1, ndim do dxm[d] = grid:dx(d) end
	    grid:cellCenter(xcm)

	    grid:setIndex(idxp)
	    for d = 1, ndim do dxp[d] = grid:dx(d) end
	    grid:cellCenter(xcp)

	    qIn:fill(qInIdxr(idxm), qInM)
	    qIn:fill(qInIdxr(idxp), qInP)

	    qOut:fill(qOutIdxr(idxm), qOutM)
	    qOut:fill(qOutIdxr(idxp), qOutP)

	    if firstDir and i<=dirUpIdx-1 then
	       local cflFreq = self._equation:volTerm(xcp, dxp, idxp, qInP, qOutP)
	       cfla = math.max(cfla, cflFreq*dt)
	    end
	    if i >= dirLoSurfIdx and i <= dirUpSurfIdx then
	       local maxs = self._equation:surfTerm(
		  dir, cfla, xcm, xcp, dxm, dxp, self._maxsOld[dir], idxm, idxp, qInM, qInP, qOutM, qOutP)
	       self._maxsLocal[dir] = math.max(self._maxsLocal[dir], maxs)
            else
	       if self._zeroFluxFlags[dir] then
	          -- we need to give equations a chance to apply partial
	          -- surface updates even when the zeroFlux BCs have been
	          -- applied
	          self._equation:boundarySurfTerm(
		     dir, xcm, xcp, dxm, dxp, self._maxsOld[dir], idxm, idxp, qInM, qInP, qOutM, qOutP)
               end
	    end
	 end
      end
      firstDir = false
   end

   -- determine largest amax across processors
  Mpi.Allreduce(
     self._maxsLocal:data(), self._maxs:data(), ndim, Mpi.DOUBLE, Mpi.MAX, self:getComm())

   -- return failure if time-step was too large
   if cfla > cflm then return false, dt*cfl/cfla end

   -- accumulate full solution if not computing increments
   if not self._onlyIncrement then
      qOut:scale(dt); qOut:accumulate(1.0, qIn) -- qOut = qIn + dt*qOut
   end

   self._isFirst = false
   return true, dt*cfl/cfla
end

return HyperDisCont
