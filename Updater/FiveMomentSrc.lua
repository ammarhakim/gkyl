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

-- Explicit, SSP RK3 scheme
local function updateSrcRk3(self, dt, fPtr, emPtr)
end

-- Use an explicit scheme to update momentum and electric field: this
-- is an extension of the standard Boris push algorithm, in which the
-- half time-step electric field update is replaced by an implicit
-- step in which both the velocity and electric are updater.
local function updateSrcModBoris(self, dt, fPtr, emPtr)
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCentered(self, dt, dt1, dt2, fPtr, emPtr)
   print("Implicit update")
end


-- Five-moment source updater object
local FiveMomentSrc = {}
-- constructor
function FiveMomentSrc:new(tbl)
   local self = setmetatable({}, FiveMomentSrc)
   Base.setup(self, tbl) -- setup base object

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

   -- scheme is one of "ssp-rk3", "modified-boris"  or "time-centered"
   local scheme = tbl.scheme and tbl.scheme or "time-centered"
   self._updateSrc = nil
   if scheme == "ssp-rk3" then
      self._updateSrc = updateSrcRk3
   elseif scheme == "modified-boris" then
      self._updateSrc = updateSrcModBoris
   elseif scheme == "time-centered" then
      self._updateSrc = updateSrcTimeCentered
   end

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

   -- make list of indexer functions and pointers
   local fIdxr, fPtr = {}, {}
   for i = 1, self._nFluids do
      fIdxr[i] = outFld[i]:genIndexer()
      fPtr[i] = outFld[i]:get(1)
   end
   local emIdxr = emFld:genIndexer()
   local emPtr = emFld:get(1)

   -- loop over local range, updating source in each cell
   local localRange = emFld:localRange()
   for idx in localRange:colMajorIter() do
      -- set pointers to fluids and field
      for i = 1, self._nFluids do
	 outFld[i]:fill(fIdxr[i](idx), fPtr[i])
      end
      emFld:fill(emIdxr(idx), emPtr)

      -- update momentum and electric field
      self._updateSrc(self, dt, fPtr, emPtr)
      -- check if we have pressure and update pressure
      
   end
end

-- Methods for wave-propagation scheme
FiveMomentSrc.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   FiveMomentSrc = FiveMomentSrc
}
