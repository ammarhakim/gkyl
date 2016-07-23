-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update five-moment source terms. This updater allows
-- both explicit and implicit updates.
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
  double _charge, _mass; /* Charge and mass */
} FluidData_t;

typedef struct {
  uint8_t _nFluids; /* Number of fluids */
  double _epsilon0; /* Permittivity of free space */
  double _chi_e, _chi_m; /* Propagation speed factor for electric field error potential */
  uint8_t _gravityDir; /* Direction of gravity force */
  double _gravity; /* Gravitational acceleration */
  uint8_t _hasStatic, _hasPressure; /* Flag to indicate if there is: static EB field, pressure */
} FiveMomentSrcData_t;
]]

-- Five-moment source updater object
local FiveMomentSrc = {}

-- constructor
function FiveMomentSrc:new(tbl)
   local self = setmetatable({}, FiveMomentSrc)

   self._onGrid = assert(tbl.onGrid, "Updater.FiveMomentSrc: Must provide grid object using 'onGrid'")   
   -- read data from input table
   self._nFluids = assert(tbl.numFluids, "Updater.FiveMomentSrc: Must specify number of fluid using 'numFluids'")

   assert(#tbl.charge == 2, "Updater.FiveMomentSrc: Charge table must have " .. self._nFluids .. " elements.")
   assert(#tbl.mass == 2, "Updater.FiveMomentSrc: Mass table must have " .. self._nFluids .. " elements.")
   self._charge, self._mass = tbl.charge, tbl.mass

   self._epsilon0 = assert(tbl.epsilon0, "Updater.FiveMomentSrc: Must specify 'epsilon0'")
   self._chi_e = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self._chi_m = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0
   
   self._hasStatic = tbl.hasStaticField and tbl.hasStaticField or false
   self._hasPressure = tbl.hasPressure and tbl.hasPressure or false
   
   self._gravity = tbl.gravity and tbl.gravity or 0.0
   self._gravityDir = tbl.dir and tbl.dir or 1 -- by default gravity acts in X direction

   -- scheme is one of "explicit" or "implicit"
   self._scheme = self.scheme and self.scheme or "implicit"

   self._qbym = new("double[?]", self._nFluids+1) -- first value is dummy to allow index from 1
   for i = 1, self._nFluids do
      self._qbym[i] = self._charge[i]/self._mass[i]
   end
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(FiveMomentSrc, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   -- check if correct number of inputs were provided
   assert(#outFld == self._nFluids+1,
	  "Must have exactly " .. self._nFluids+1 .. " output fields. Provided " .. #inFld .. " instead.")
   local emFld = outFld[#outFld] -- EM field

   local dt1 = 0.5*dt
   local dt2 = 0.5*dt/self._epsilon0

   -- loop over local range, updating source in each cell
   local localRange = emFld:localRange()
   for idx in localRange:colMajorIter() do
      
   end
end

-- Methods for wave-propagation scheme
FiveMomentSrc.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   FiveMomentSrc = FiveMomentSrc
}
