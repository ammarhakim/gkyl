-- Gkyl ------------------------------------------------------------------------
--
-- Finite volume wave propagation updater on rectangular, Cartesian
-- grid. Non-uniform grids are supported and an embedded BC can also
-- be included. Algorithm based on LeVeque book and Hakim
-- et. al. papers. Original design from my (AHH) thesis code,
-- Miniwarpx and Warpx which significantly extend ideas present in
-- CLAWPACK. This implementation in LuaJIT is, perhaps, completely
-- novel and unique.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
"new, copy, fill, sizeof, typeof, metatype")

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
   waves1[${MEQN}*mw1+${i}]*waves[${MEQN}*mw1+${i}]+
   |end
   waves1[${MEQN}*mw1+${MEQN-1}]*waves[${MEQN}*mw1+${MEQN-1}]
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

-- limiter functions: See LeVeque book
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
   local c = (1.0+r)/2
   return math.max(0.0, math.min(c, 2, 2*r))
end
limiterFunctions["beam-warming"] = function (r)
   return r
end
limiterFunctions["zero"] = function (r)
   return 0
end

-- Helper object for indexing 1D slice data. The slice spans from
-- [lower, upper] (inclusive) and has `stride` pieces of data stored
-- at each location.
local createSliceData = function (dtype)
   local slice_mt = {
      __new = function (self, lower, upper, stride)
         local n = upper-lower+1
         local v = new(self)
         v._data = Alloc.malloc(typeof(dtype)*n*stride)
         v._sz, v._stride, v._lower = n*stride, stride, lower
         return v
      end,
      __index = function (self, k)
         return self._data+(k-self._lower)*self._stride
      end,
      __gc = function (self)
         Alloc.free(self._data)
      end,
   }
   return metatype(typeof(string.format("struct {int32_t _sz, _stride, _lower; %s *_data; }", dtype)), slice_mt)
end
local SliceData = createSliceData("double")
local SliceDataBool = createSliceData("bool")

-- helper function to zero out contents of SliceData
local function clearSliceData(sd, dtype)
   fill(sd._data, sd._sz*sizeof("double"))
end

-- Wave-propagation updater object
local WavePropagation = Proto(UpdaterBase)

-- constructor
function WavePropagation:init(tbl)
   WavePropagation.super.init(self, tbl) -- setup base object

   self._isFirst = true -- will be reset first time _advance() is called

   -- read data from input table
   self._onGrid = assert(tbl.onGrid, "Updater.WavePropagation: Must provide grid object using 'onGrid'")
   self._equation = assert(tbl.equation, "Updater.WavePropagation: Must provide equation object using 'equation'")
   self._limiterFunc = assert(limiterFunctions[tbl.limiter], "Updater.WavePropagation: Must specify limiter to use")
   self._limiterFuncZero = limiterFunctions["zero"]

   self._ndim = self._onGrid:ndim()
   self._cfl = assert(tbl.cfl, "Updater.WavePropagation: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl
   self._hasSsBnd = xsys.pickBool(tbl.hasSsBnd, false)
   self._inOut = tbl.inOut

   self._updateDirs = {}
   local upDirs = tbl.updateDirections and tbl.updateDirections or {1, 2, 3, 4, 5, 6}
   for d = 1, #upDirs do
      self._updateDirs[d] = upDirs[d] -- update directions
   end

   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()
   self._localRange = tbl.onGrid:localRange()
   local localRange = self._localRange

   -- allocate space for storing 1D slice data
   self.wavesSlice, self.speedsSlice, self.fsSlice = {}, {}, {}
   for d = 1, self._ndim do
      local l, u = localRange:lower(d)-2, localRange:upper(d)+2
      self.wavesSlice[d] = SliceData(l, u, meqn*mwave)
      self.speedsSlice[d] = SliceData(l, u, mwave)
      self.fsSlice[d] = SliceData(l, u, meqn)
   end

   -- masks for edge/cell locations regarding stari-stepped boundaries
   self.onSsBnd = {}  -- if exactly one of the edge's left/right cells is outside
   self.bothOutSsBnd = {} -- if both of the edge's left/right cells are outside
   self.thisOutSsBnd = {}  -- if this cell is outside
   if (self._hasSsBnd) then
      for d = 1, self._ndim do
         local l, u = localRange:lower(d)-2, localRange:upper(d)+2
         self.onSsBnd[d] = SliceDataBool(l, u, 1)
         self.bothOutSsBnd[d] = SliceDataBool(l, u, 1)
         self.thisOutSsBnd[d] = SliceDataBool(l, u, 1)
      end
   end

   -- store range objects needed in update
   self._perpRangeDecomp = {}

   -- construct various functions from template representations
   self._calcDelta = loadstring( calcDeltaTempl {MEQN = meqn} )()
   self._calcCfla = loadstring( calcCflaTempl {MWAVE = mwave} )()
   self._calcFirstOrderGud = loadstring( calcFirstOrderGudTempl {MEQN = meqn} )()
   self._waveDotProd = loadstring( waveDotProdTempl {MEQN = meqn} )()
   self._rescaleWave = loadstring( rescaleWaveTempl {MEQN = meqn} )()
   self._secondOrderFlux = loadstring( secondOrderFluxTempl {MEQN = meqn} )()
   self._secondOrderUpdate = loadstring( secondOrderUpdateTempl {MEQN = meqn} )()
end

-- Limit waves: this code closely follows the example of CLAWPACK and
-- my (AHH) thesis code Miniwarpx.
function WavePropagation:limitWaves(lower, upper, wavesSlice, speedsSlice, onSsBnd)
   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()
   for mw = 1, mwave do
      local dotr = self._waveDotProd(meqn, wavesSlice[lower-1], wavesSlice[lower], mw)
      for i = lower, upper do
         local dotl = dotr
         local wnorm2 = self._waveDotProd(meqn, wavesSlice[i], wavesSlice[i], mw)
         dotr = self._waveDotProd(meqn, wavesSlice[i], wavesSlice[i+1], mw)
         if wnorm2 > 0 then
            local r = speedsSlice[i][mw-1] > 0 and dotl/wnorm2 or dotr/wnorm2
            local wlimitr
            if (self._hasSsBnd and onSsBnd[i][0]) then
               wlimitr = self._limiterFuncZero(r)
            else
               wlimitr = self._limiterFunc(r)
            end
            self._rescaleWave(wlimitr, wavesSlice[i]+(mw-1)*meqn)
         end
      end
   end
end

local isOutside = function (inOutPtr)
   return inOutPtr[0] < 0.
end

-- advance method
function WavePropagation:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local qIn = assert(inFld[1], "WavePropagation.advance: Must-specify an input field")
   local qOut = assert(outFld[1], "WavePropagation.advance: Must-specify an output field")

   local equation = self._equation -- equation to solve
   local meqn, mwave = equation:numEquations(), equation:numWaves()
   local localRange = self._localRange

   local qInIdxr, qOutIdxr = qIn:genIndexer(), qOut:genIndexer() -- indexer functions into fields

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local delta = Lin.Vec(meqn)
   local waves = Lin.Mat(mwave, meqn)
   local s = Lin.Vec(mwave)
   local amdq, apdq = Lin.Vec(meqn), Lin.Vec(meqn)

   local idxp, idxm = Lin.IntVec(grid:ndim()), Lin.IntVec(grid:ndim())

   -- pointers to (re)use inside inner loop
   local qInL, qInR = qIn:get(1), qIn:get(1)
   local qOutL, qOutR = qOut:get(1), qOut:get(1)
   local q1 = qOut:get(1)

   local inOutL, inOutR
   local inOut1
   local inOutIdxr
   if (self._hasSsBnd) then
      inOutL, inOutR = self._inOut:get(1), self._inOut:get(1)
      inOut1 = self._inOut:get(1)
      inOutIdxr = self._inOut:genIndexer()
   end

   local tId = grid:subGridSharedId() -- local thread ID

   qOut:copy(qIn) -- update only adds increments, so set qOut = qIn
   -- update specified directions
   for d = 1, #self._updateDirs do
      local dir = self._updateDirs[d]
      local dtdx = dt/grid:dx(dir)

      local wavesSlice, speedsSlice, fsSlice = self.wavesSlice[dir], self.speedsSlice[dir], self.fsSlice[dir]
      local onSsBnd = self.onSsBnd[dir]
      local bothOutSsBnd = self.bothOutSsBnd[dir]
      local thisOutSsBnd = self.thisOutSsBnd[dir]
      -- lower/upper bounds in direction 'dir': these are edge indices
      local dirLoIdx, dirUpIdx = localRange:lower(dir)-1, localRange:upper(dir)+2

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
      for idx in perpRangeDecomp:colMajorIter(tId) do
         idx:copyInto(idxp); idx:copyInto(idxm)

         -- fill masks along this direction only once
         if self._isFirst and self._hasSsBnd then
            for i = localRange:lower(dir)-1, localRange:upper(dir)+2 do -- this loop is over edges
               idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'
               self._inOut:fill(inOutIdxr(idxm), inOutL)
               self._inOut:fill(inOutIdxr(idxp), inOutR)
               local isOutsideL = isOutside(inOutL)
               local isOutsideR = isOutside(inOutR)
               onSsBnd[i][0] =  (isOutsideL and not isOutsideR) or
                  (not isOutsideL and isOutsideR)
               bothOutSsBnd[i][0] =  isOutsideL and isOutsideR
               thisOutSsBnd[i][0] = isOutsideL
            end
         end

         for i = dirLoIdx, dirUpIdx do -- this loop is over edges
            idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'

            if not (self._hasSsBnd and bothOutSsBnd[i][0]) then
               qIn:fill(qInIdxr(idxm), qInL); qIn:fill(qInIdxr(idxp), qInR)
               self._calcDelta(qInL, qInR, delta) -- jump across interface

               equation:rp(dir, delta, qInL, qInR, waves, s) -- compute waves and speeds
               equation:qFluctuations(dir, qInL, qInR, waves, s, amdq, apdq) -- compute fluctuations

               qOut:fill(qOutIdxr(idxm), qOutL); qOut:fill(qOutIdxr(idxp), qOutR)
               self._calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq) -- first-order Gudonov updates
               cfla = self._calcCfla(cfla, dtdx, s) -- actual CFL value

               -- copy waves data for use in limiters
               copy(wavesSlice[i], waves:data(), sizeof("double")*meqn*mwave)
               copy(speedsSlice[i], s:data(), sizeof("double")*mwave)
            end
         end

         -- return if time-step was too large
         if cfla > cflm then return false, dt*cfl/cfla end

         -- limit waves before computing second-order updates
         self:limitWaves(localRange:lower(dir), localRange:upper(dir)+1, wavesSlice, speedsSlice, onSsBnd)

         local dirLoIdx2, dirUpIdx2 = localRange:lower(dir), localRange:upper(dir)+1 -- one more edge than cells
         -- compute second order correction fluxes
         clearSliceData(fsSlice)
         for i = dirLoIdx2, dirUpIdx2 do -- this loop is over edges
            if not (self._hasSsBnd and bothOutSsBnd[i][0]) then
               for mw = 0, mwave-1 do
                  self._secondOrderFlux(dtdx, speedsSlice[i][mw], wavesSlice[i]+meqn*mw, fsSlice[i])
               end
            end
         end

         local dirLoIdx3, dirUpIdx3 = localRange:lower(dir), localRange:upper(dir)
         -- add them to solution
         for i = dirLoIdx3, dirUpIdx3 do -- this loop is over cells
            idxm[dir] = i -- cell index
            if not (self._hasSsBnd and thisOutSsBnd[i][0]) then
               qOut:fill(qOutIdxr(idxm), q1)
               self._secondOrderUpdate(dtdx, fsSlice[i], fsSlice[i+1], q1)
            end
         end
      end
   end

   self._isFirst = false -- no longer first time
   return true, dt*cfl/cfla
end

return WavePropagation
