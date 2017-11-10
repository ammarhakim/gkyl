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
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"
local _ = require "Updater.linearHyperbolicData.LinearHyperModDecl"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- pick proper function for volume integrals
local function selectVolFunc(bId, polyOrder, ndim, dir)
   if bId == "serendipity" then
      return ffi.C[string.format("LinearHyperSer%dDP%d_Vol%d", ndim, polyOrder, dir)]
   elseif bId == "maximal-order" then
      return ffi.C[string.format("LinearHyperMax%dDP%d_Vol%d", ndim, polyOrder, dir)]
   end
end

-- Linear hyperbolic DG solver updater object
local LinearHyperbolicDisCont = {}

function LinearHyperbolicDisCont:new(tbl)
   local self = setmetatable({}, LinearHyperbolicDisCont)
   Base.setup(self, tbl) -- setup base object

   -- read data from input file
   self._onGrid = assert(tbl.onGrid, "Updater.LinearHyperbolicDisCont: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.LinearHyperbolicDisCont: Must specify basis functions to use using 'basis'")

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   self._ndim = self._onGrid:ndim()

   -- linear equation to solve
   self._equation = assert(tbl.equation, "Updater.LinearHyperbolicDisCont: Must provide equation object using 'equation'")

   -- read in which directions we are to update
   self._nUpdateDirs = tbl.updateDirections and #tbl.updateDirections or self._ndim
   local upDirs = tbl.updateDirections and tbl.updateDirections or {1, 2, 3}
   self._updateDirs = {}
   for d = 1, #upDirs do
      self._updateDirs[d] = upDirs[d] -- update directions
   end   

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.LinearHyperbolicDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- functions to perform volume updates
   self._volumeUpdate = {}
   for d = 1, self._ndim do
      self._volumeUpdate[d] = selectVolFunc(self._basis:id(), self._basis:polyOrder(), self._ndim, d)
   end
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(LinearHyperbolicDisCont, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   local grid = self._onGrid

   local qIn = assert(outFld[1], "LinearHyperbolicDisCont.advance: Must specify an input field")
   local qOut = assert(outFld[1], "LinearHyperbolicDisCont.advance: Must specify an output field")

   local ndim = grid:ndim()
   local numBasis = self._basis:numBasis()
   local numSurfBasis = self._basis:numSurfBasis()
   local meqn = qOut:numComponents()/numBasis

   local flux = Lin.Vec(qOut:numComponents()) -- flux 
   local numericalFlux = Lin.Vec(numSurfBasis*meqn) -- numerical flux on face 

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local localRange = qOut:localRange()
   local qInIdxr, qOutIdxr = qIn:genIndexer(), qOut:genIndexer() -- indexer functions into fields

   -- to store grid info
   local dx = Lin.Vec(ndim) -- cell shape
   local xc = Lin.Vec(ndim) -- cell center

   qOut:clear(0.0) -- compute increments initally

   -- pointers for (re)use in update
   local qInPtr, qOutPtr = qIn:get(1), qOut:get(1)
   local qInL, qInR = qIn:get(1), qIn:get(1)
   local qOutL, qOutR = qOut:get(1), qOut:get(1)   
   
   -- accumulate contributions from volume integrals
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      for d = 1, ndim do dx[d] = grid:dx(d) end
      local volFact = grid:cellVolume()/2^ndim

      qIn:fill(qInIdxr(idx), qInPtr)
      qOut:fill(qOutIdxr(idx), qOutPtr)
      
      for d = 1, self._nUpdateDirs do -- update only specified directions
	 local dir = self._updateDirs[d]
	 self._equation:fluxCoeff(dir, basis, qInPtr, flux) -- compute flux
	 self._volumeUpdate[dir](meqn, numBasis, volFact*2/dx[dir], flux:data(), qOutPtr:data())
      end
   end

   -- accumulate contributions from surface integrals
   for d = 1, self._nUpdateDirs do
      local dir = self._updateDirs[d]

      -- lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell)
      local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`.
      for idx in perpRange:colMajorIter() do
	 local idxp, idxm = idx:copy(), idx:copy()

   	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'

	    qIn:fill(qInIdxr(idxm), qInL)
	    qIn:fill(qInIdxr(idxp), qInR)

	    -- compute left/right expansions on faces

	    -- compute numerical flux at the common face
	    
	    -- update two cells connected to face

	    -- compute actual CFL number used

	 end
      end
   end
   
   return true, GKYL_MAX_DOUBLE
end

-- Methods in updater
LinearHyperbolicDisCont.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   LinearHyperbolicDisCont = LinearHyperbolicDisCont
}
