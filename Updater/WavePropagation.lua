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

-- Template for function to compute first-order Gudonov update
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

-- Template for function to compute dot product of waves
local waveDotProdTempl = xsys.template([[
return function (meqn, waves, waves1, mw)
  local mw1 = mw-1
  return
|for i = 0, MEQN-2 do
  waves[meqn*mw1+${i}]*waves[meqn*mw1+${i}]+
|end
  waves1[meqn*mw1+${MEQN-1}]*waves[meqn*mw1+${MEQN-1}]
end
]])

-- Template for function to rescale waves
local rescaleWaveTempl = xsys.template([[
return function (scale, wave)
|for i = 0, MEQN-1 do
  wave[${i}] = scale*wave[${i}]
|end
end
]])

local secondOrderFluxTempl = xsys.template([[
return function (dtdx, s, wave, fs)
  local sfact = 0.5*math.abs(s)*(1-math.abs(s)*dtdx)
|for i = 0, MEQN-1 do
  fs[${i}] = fs[${i}] + sfact*wave[${i}]
|end
end
]])

local secondOrderUpdateTempl = xsys.template([[
return function (dtdx, fs, fs1, q)
| for i = 1, MEQN do
  q[${i}] = q[${i}] - dtdx*(fs1[${i-1}]-fs[${i-1}])
| end
end
]])

-- limiter functions
local limiterFunctions = {}
limiterFunctions["no-limiter"] = function (r)
   return 1
end
limiterFunctions["min-mod"] = function (r)
   return math.max(0, math.min(1, r))
end
limiterFunctions["superbee"] = function (r)
   return math.max(0.0, math.min(1, 2*r), math.min(2.0, r))
end
limiterFunctions["van-leer"] = function (r)
   return (r+math.abs(r))/(1+math.abs(r))
end
limiterFunctions["monotonized-centered"] = function (r)
   return (1.0+r)/2
end
limiterFunctions["beam-warming"] = function (r)
   return r
end
limiterFunctions["beam-warming"] = function (r)
   return 0
end

-- Helper object for indexing 1D slice data. This is zero-indexed as
-- the underlying data is just a raw-C array
local slice_mt = {
   __new = function (self, nsize, stride)
      local v = new(self, nsize*stride)
      v._sz, v._stride = nsize*stride, stride
      return v
   end,
   __index = function (self, k)
      return self._data+k*self._stride
   end,
}
local SliceData = metatype(typeof("struct {int32_t _sz, _stride; double _data[?]; }"), slice_mt)
-- helper function to zero out contents of SliceData
function clearSliceData(sd)
   fill(sd._data, sd._sz*sizeof("double"))
end

-- Wave-propagation updater object
WavePropagation = {}

-- constructor
function WavePropagation:new(tbl)
   local self = setmetatable({}, WavePropagation)

   -- read data from input table
   self._onGrid = assert(tbl.onGrid, "Updater.WavePropagation: Must provide grid object using 'onGrid'")
   self._equation = assert(tbl.equation, "Updater.WavePropagation: Must provide equation object using 'equation'")
   self._limiterFunc = assert(limiterFunctions[tbl.limiter], "Updater.WavePropagation: Must specify limiter to use")

   -- set private data
   self._privData = new(typeof("WavePrivateData_t"))
   self._privData._ndim = self._onGrid:ndim()
   self._privData._cfl = assert(tbl.cfl, "Updater.WavePropagation: Must specify CFL number using 'cfl'")
   self._privData._cflm = tbl.cflm and tbl.cflm or 1.1*self._privData._cfl

   self._privData._nUpdateDirs = tbl.updateDirections and #tbl.updateDirections or self._privData._ndim
   local upDirs = tbl.updateDirections and tbl.updateDirections or {1, 2, 3, 4, 5, 6}
   for d = 1, self._privData._ndim do
      self._privData._updateDirs[d] = upDirs[d] -- update directions
   end

   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()

   -- allocate space for storing 1D slice data   
   local shapeMax = 0
   local localRange = tbl.onGrid:localRange()
   for d = 1, self._privData._ndim do
      shapeMax = math.max(shapeMax, localRange:shape(d))
   end
   local sz = (shapeMax+4) -- 2 ghost cells on each side
   self.wavesSlice = SliceData(sz, tbl.equation:numWaves()*tbl.equation:numEquations())
   self.speedsSlice = SliceData(sz, tbl.equation:numWaves())
   self.fsSlice = SliceData(sz, tbl.equation:numEquations())

   -- construct various functions from template representations
   self._calcDelta = loadstring( calcDeltaTempl {MEQN = meqn} )()
   self._calcCfla = loadstring( calcCflaTempl {MWAVE = mwave} )()
   self._calcFirstOrderGud = loadstring( calcFirstOrderGudTempl {MEQN = meqn} )()
   self._waveDotProd = loadstring( waveDotProdTempl {MEQN = meqn} )()
   self._rescaleWave = loadstring( rescaleWaveTempl {MEQN = meqn} )()
   self._secondOrderFlux = loadstring( secondOrderFluxTempl {MEQN = meqn} )()
   self._secondOrderUpdate = loadstring( secondOrderUpdateTempl {MEQN = meqn} )()

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(WavePropagation, { __call = function (self, o) return self.new(self, o) end })

-- Limit waves: this code closely follows the example of CLAWPACK and
-- my (AHH) thesis code Miniwarpx.
local function limitWaves(self, limitRange, wavesSlice, speedsSlice)
   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()
   for mw = 1, mwave do
      local dotr = self._waveDotProd(meqn, wavesSlice[0], wavesSlice[1], mw)
      for i = 1, limitRange do
	 local dotl = dotr
	 local wnorm2 = self._waveDotProd(meqn, wavesSlice[i], wavesSlice[i], mw)
	 dotr = self._waveDotProd(meqn, wavesSlice[i], wavesSlice[i+1], mw)
	 if wnorm2 > 0 then
	    local r = speedsSlice[i][mw-1] > 0 and dotl/wnorm2 or dotr/wnorm2
	    local wlimitr = self._limiterFunc(r)
	    self._rescaleWave(wlimitr, wavesSlice[i])
	 end
      end
   end
end

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
   local wavesStride, speedsStride = mwave*meqn, mwave

   -- update specified directions
   for d = 1, self._privData._nUpdateDirs do
      local dir = self._privData._updateDirs[d]
      local dtdx = dt/grid:dx(dir)
      
      -- lower/upper bounds in direction 'dir': these are edge indices
      local dirLoIdx, dirUpIdx = localRange:lower(dir)-1, localRange:upper(dir)+2
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'
      -- for wave limiters
      local limitRange = localRange:shape(dir)+1 -- there is one more edge than cells

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`. 
      for idx in perpRange:colMajorIter() do
	 local idxp, idxm = idx:copy(), idx:copy()

	 local count = 0 --location into the 1D slice storage
   	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'
	    
	    local qInL, qInR = qIn:get(qInIdxr(idxm)), qIn:get(qInIdxr(idxp))
	    self._calcDelta(qInL, qInR, delta) -- jump across interface

	    equation:rp(delta, qInL, qInR, waves, s) -- compute waves and speeds from jump
	    equation:qFluctuations(qInL, qInR, waves, s, amdq, apdq) -- compute fluctuations

	    local qOutL, qOutR = qOut:get(qOutIdxr(idxm)), qOut:get(qOutIdxr(idxp))
	    self._calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq) -- first-order Gudonov updates

	    cfla = self._calcCfla(cfla, dtdx, s) -- actual CFL value

	    -- copy waves data for use in limiters
	    copy(self.wavesSlice[count], waves:data(), sizeof("double")*meqn*mwave)
	    copy(self.speedsSlice[count], s:data(), sizeof("double")*mwave)

	    count = count+1
	 end
	 -- return if time-step was too large
	 if cfla > cflm then return false, dt*cfl/cfla end

	 -- limit waves before computing second-order updates
	 limitWaves(self, limitRange, self.wavesSlice, self.speedsSlice)

	 local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1 -- one more edge than cells
	 clearSliceData(self.fsSlice)
	 -- compute second order correction fluxes
	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    for mw = 0, mwave-1 do
	       self._secondOrderFlux(dtdx, self.speedsSlice[i][mw], self.wavesSlice[i]+meqn*mw, self.fsSlice[i])
	    end
	 end

	 local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)
	 -- add them to solution
	 for i = dirLoIdx, dirUpIdx do -- this loop is over cells
	    idxm[dir] = i -- cell index
	    local q1 = qOut:get(qOutIdxr(idxm))
	    --self._secondOrderUpdate(dtdx, self.fsSlice[i], self.fsSlice[i+1], q1)
	 end
      end
   end
   return true, dt*cfl/cfla
end

-- Methods for wave-propagation scheme methods
WavePropagation.__index = {
   advance = Base.advanceFuncWrap(advance)
}

return {
   WavePropagation = WavePropagation
}
