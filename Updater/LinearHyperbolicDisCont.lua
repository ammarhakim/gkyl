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

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Linear hyperbolic DG solver updater object
local LinearHyperbolicDisCont = {}

function LinearHyperbolicDisCont:new(tbl)
   local self = setmetatable({}, LinearHyperbolicDisCont)
   Base.setup(self, tbl) -- setup base object

   -- read data from input file
   self._onGrid = assert(tbl.onGrid, "Updater.LinearHyperbolicDisCont: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.LinearHyperbolicDisCont: Must specify basis functions to use using 'basis'")

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

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
   local meqn = qOut:numComponents()/numBasis

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local localRange = qOut:localRange()
   local indexer = qOut:genIndexer()

   -- to store grid info
   local dx = Lin.Vec(ndim) -- cell shape
   local xc = Lin.Vec(ndim) -- cell center   

   -- pointers for (re)use in update
   local qInPtr, qOutPtr = qIn:get(1), qOut:get(1)
   local qInL, qInR = qIn:get(1), qIn:get(1)
   local qOutL, qOutR = qOut:get(1), qOut:get(1)
   
   -- Accumulate contributions from volume integrals
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      -- get cell shape, cell center coordinates
      for d = 1, ndim do dx[d] = grid:dx(d) end

      -- set pointers
      q
      
      -- update only specified directions
      for d = 1, self._nUpdateDirs do
	 local dir = self._updateDirs[d]

	 
      end
      
   end
   
   return true, GKYL_MAX_DOUBLE
end

-- Methods in updater
LinearHyperbolicDisCont.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   LinearHyperbolicDisCont = LinearHyperbolicDisCont
}
