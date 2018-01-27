-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute RHS or forward Euler update for linear
-- hyperbolic equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local ffi = require "ffi"
local xsys = require "xsys"

-- Hyperbolic DG solver updater object
local HyperDisCont = Proto(UpdaterBase)

function HyperDisCont:init(tbl)
   HyperDisCont.super.init(self, tbl)

   self._isFirst = true -- will be reset first time _advance() is called

   -- read data from input file
   self._onGrid = assert(tbl.onGrid, "Updater.HyperDisCont: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.HyperDisCont: Must specify basis functions to use using 'basis'")

   -- by default, compute forward Euler: if onlyIncrement is true,
   -- then only increments are computed. NOTE: The increments are NOT
   -- multiplied by dt.
   self._onlyIncrement = xsys.pickBool(tbl.onlyIncrement, false)

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

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.HyperDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- store range objects needed in update
   self._perpRange = {}

   return self
end

-- advance method
function HyperDisCont:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid

   local qIn = assert(inFld[1], "HyperDisCont.advance: Must specify an input field")
   local qOut = assert(outFld[1], "HyperDisCont.advance: Must specify an output field")

   local ndim = grid:ndim()

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local localRange = qOut:localRange()
   local qInIdxr, qOutIdxr = qIn:genIndexer(), qOut:genIndexer() -- indexer functions into fields

   -- to store grid info
   local dx = Lin.Vec(ndim) -- cell shape
   local xc = Lin.Vec(ndim) -- cell center
   local idxp, idxm = Lin.IntVec(ndim), Lin.IntVec(ndim)

   -- pointers for (re)use in update
   local qInM, qInP = qIn:get(1), qIn:get(1)
   local qOutM, qOutP = qOut:get(1), qOut:get(1)

   -- This flag is needed as the volume integral already contains
   -- contributions from all directions. Hence, we can only
   -- accumulate the volume contribution once, skipping it for
   -- other directions
   local firstDir = true

   qOut:clear(0.0) -- compute increments
   -- accumulate contributions from surface integrals
   for _, dir in ipairs(self._updateDirs) do
      -- lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell)
      local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1

      if self._isFirst then
	 self._perpRange[dir] = localRange:shorten(dir) -- range orthogonal to 'dir'
      end
      local perpRange = self._perpRange[dir]

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`.
      for idx in perpRange:colMajorIter() do
	 idx:copyInto(idxp); idx:copyInto(idxm)
   	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'
	    -- compute cell center coordinates and cell spacing
	    grid:setIndex(idxp)
	    for d = 1, ndim do dx[d] = grid:dx(d) end
	    grid:cellCenter(xc)

	    qIn:fill(qInIdxr(idxm), qInM)
	    qIn:fill(qInIdxr(idxp), qInP)

	    qOut:fill(qOutIdxr(idxm), qOutM)
	    qOut:fill(qOutIdxr(idxp), qOutP)

	    local maxs = self._equation:maxSpeed(dir, xc, dx, qInP)
	    cfla = math.max(cfla, maxs*dt/dx[dir])

	    if firstDir then
	       self._equation:volTerm(xc, dx, qInP, qOutP)
	    end
	    self._equation:surfTerm(dir, xc, dx, qInM, qInP, qOutM, qOutP)
	 end
	 -- return failure if time-step was too large
	 if cfla > cflm then return false, dt*cfl/cfla end
      end
      firstDir = false
   end

   -- accumulate full solution if not computing increments
   if not self._onlyIncrement then
      qOut:scale(dt); qOut:accumulate(1.0, qIn) -- qOut = qIn + dt*qOut
   end

   self._isFirst = false
   return true, dt*cfl/cfla
end

return HyperDisCont
