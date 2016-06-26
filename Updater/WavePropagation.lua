-- Gkyl ------------------------------------------------------------------------
--
-- Finite volume wave propagation updater. Algorithm based on LeVeque
-- book and Hakim et. al. papers. Original design from my (AHH) thesis
-- code, Miniwarpx and Warpx which significantly extend ideas present
-- in CLAWPACK. This implementation in LuaJIT is, perhaps, completely
-- novel and unique.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")


-- Define C types for storing private data for use in updater
ffi.cdef [[
typedef struct {
    uint8_t _ndim; /* Solver dimension */
    double _cfl, _cflm; /* CFL number and maximum-CFL number */
    uint8_t _nUpdateDirs, _updateDirs[6]; /* Number and directions to update */
} WavePrivateData_t ;
]]

-- Template for function to compute jump 
local calcDeltaTempl = xsys.template([[
return function (ql, qr, delta)
|for i = 1, MEQN do
  delta[${i}] = qr[${i}] - ql[${i}]
|end
end
]])

-- Template for function to compute maximum CFL number
local calcCflaTempl = xsys.template([[
return function (cfla, dtdx, s)
  local c = cfla
|for i = 1, MWAVE do
  c = math.max(c, dtdx*math.abs(s[${i}]))
|end
  return c
end
]])

local calcFirstOrderGudTempl = xsys.template([[
return function (dtdx, ql, qr, amdq, apdq)
|for i = 1, MEQN do
  qr[${i}] = qr[${i}] - dtdx*apdq[${i}]
|end
|for i = 1, MEQN do
  ql[${i}] = ql[${i}] - dtdx*amdq[${i}]
|end
end
]])

-- Wave-propagation updater object
WavePropagation = {}

-- constructor
function WavePropagation:new(tbl)
   local self = setmetatable({}, WavePropagation)

   -- read data from input table
   self._onGrid = assert(tbl.onGrid, "Updater.WavePropagation: Must provide grid object using 'onGrid'")
   self._equation = assert(tbl.equation, "Updater.WavePropagation: Must provide equation object using 'equation'")
   local limiterStr = tbl.limiter and tbl.limiter or "no-limiter"

   -- set private data
   self._privData = new(typeof("WavePrivateData_t"))
   self._privData._ndim = self._onGrid:ndim()
   self._privData._cfl = assert(tbl.cfl, "Updater.WavePropagation: Must specify CFL number using 'cfl'")
   self._privData._cflm = 1.1*self._privData._cfl

   self._privData._nUpdateDirs = tbl.updateDirections and #tbl.updateDirections or self._privData._ndim
   local upDirs = tbl.updateDirections and tbl.updateDirections or {1, 2, 3, 4, 5, 6}
   for d = 1, self._privData._ndim do
      self._privData._updateDirs[d] = upDirs[d] -- update directions
   end

   -- allocate space for storing 1D slice data
   

   -- construct various functions from template representations
   self._calcDelta = loadstring( calcDeltaTempl {MEQN = self._equation:numEquations()} )()
   self._calcCfla = loadstring( calcCflaTempl {MWAVE = self._equation:numWaves()} )()
   self._calcFirstOrderGud = loadstring( calcFirstOrderGudTempl {MEQN = self._equation:numEquations()} )()

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(WavePropagation, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local qIn = assert(inFld[1], "WavePropagation.advance: Must-specify an input field")
   local qOut = assert(outFld[1], "WavePropagation.advance: Must-specify an output field")

   local equation = self._equation -- equation to solve
   local meqn, mwave = equation:numEquations(), equation:numWaves()
   local localRange = qIn:localRange()

   local qInIdxr, qOutIdxr = qIn:genIndexer(), qOut:genIndexer() -- indexer functions into fields

   local cfl, cflm = self._privData._cfl, self._privData._cflm
   local cfla = 0.0 -- actual CFL number used
   
   local delta = Lin.Vec(meqn)
   local waves, s = Lin.Mat(mwave, meqn), Lin.Vec(mwave)
   local amdq, apdq = Lin.Vec(meqn), Lin.Vec(meqn)

   -- update specified directions
   for d = 1, self._privData._nUpdateDirs do
      local dir = self._privData._updateDirs[d]
      local dtdx = dt/grid:dx(dir)
      
      -- lower and upper bounds to loop over in direction 'dir'
      local dirLoIdx, dirUpIdx = localRange:lower(dir)-1, localRange:upper(dir)+2
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`. 
      for idx in perpRange:colMajorIter() do
	 local idxp, idxm = idx:copy(), idx:copy()

   	 for i = dirLoIdx, dirUpIdx do	    
	    idxm[dir], idxp[dir]  = i-1, i -- left/right of edge 'i'
	    
	    local qInL, qInR = qIn:get(qInIdxr(idxm)), qIn:get(qInIdxr(idxp))
	    self._calcDelta(qInL, qInR, delta) -- jump across interface

	    equation:rp(delta, qInL, qInR, waves, s) -- compute waves and speeds from jump
	    equation:qFluctuations(qInL, qInR, waves, s, amdq, apdq) -- compute fluctuations

	    local qOutL, qOutR = qOut:get(qOutIdxr(idxm)), qOut:get(qOutIdxr(idxp))
	    -- first-order Gudonov updates
	    self._calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq)

	    cfla = self._calcCfla(cfla, dtdx, s) -- actual CFL value
	 end
	 -- return if time-step was too large
	 if cfla > cflm then return false, dt*cfl/cfla end
      end
   end
   return true, dt*cfl/cfla
end

WavePropagation.__index = {
   advance = Base.advanceFuncWrap(advance)
}

return {
   WavePropagation = WavePropagation
}
