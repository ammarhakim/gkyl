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


-- Wave-propagation updater object
WavePropagation = {}

-- constructor
function WavePropagation:new(tbl)
   local self = setmetatable({}, WavePropagation)

   -- read input from input table
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
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(WavePropagation, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   -- fetch grid and input/output fields
   local grid = self._onGrid
   local qIn = assert(inFld[1], "WavePropagation.advance: Must-specify an input field")
   local qOut = assert(outFld[1], "WavePropagation.advance: Must-specify an output field")

   local equation = self._equation -- equation to solve
   local meqn, mwave = equation:numEquations(), equation:numWaves()
   local localRange = qIn:localRange()

   for d = 1, self._privData._nUpdateDirs do
      local dir = self._privData._updateDirs[d]
      
   end

   return status, 0.001
end

WavePropagation.__index = {
   advance = Base.advanceFuncWrap(advance)
}

return {
   WavePropagation = WavePropagation
}
