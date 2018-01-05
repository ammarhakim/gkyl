-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update five-moment source terms. This updater allows
-- both explicit and implicit updates.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

-- system libraries
local ffi = require "ffi"

-- Define C types for storing private data for use in updater
ffi.cdef [[
typedef struct {
  double charge, mass; /* Charge and mass */
} FluidData_t;

typedef struct {
  int8_t nFluids; /* Number of fluids */
  double epsilon0; /* Permittivity of free space */
  double chi_e, chi_m; /* Propagation speed factor for electric field error potential */
  int8_t gravityDir; /* Direction of gravity force */
  double gravity; /* Gravitational acceleration */
  bool hasStatic, hasPressure; /* Flag to indicate if there is: static EB field, pressure */
} FiveMomentSrcData_t;

  void gkylFiveMomentSrcRk3(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em);
]]

-- Explicit, SSP RK3 scheme
local function updateSrcRk3(self, dt, fPtr, emPtr)
   ffi.C.gkylFiveMomentSrcRk3(self._sd, self._fd, dt, fPtr, emPtr)
end

-- Use an explicit scheme to update momentum and electric field: this
-- is an extension of the standard Boris push algorithm, in which the
-- half time-step electric field update is replaced by an implicit
-- step in which both the velocity and electric are updated.
local function updateSrcModBoris(self, dt, fPtr, emPtr)
   print("updateSrcModBoris")
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCentered(self, dt, dt1, dt2, fPtr, emPtr)
   print("updateSrcTimeCentered")
end

-- Five-moment source updater object
local FiveMomentSrc = Proto(UpdaterBase)

function FiveMomentSrc:init(tbl)
   FiveMomentSrc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.FiveMomentSrc: Must provide grid object using 'onGrid'")

   self._sd = ffi.new(ffi.typeof("FiveMomentSrcData_t"))   
   -- Read in solver parameters
   self._sd.nFluids = assert(tbl.numFluids, "Updater.FiveMomentSrc: Must specify number of fluid using 'numFluids'")

   assert(#tbl.charge == self._sd.nFluids, "Updater.FiveMomentSrc: Charge table must have " .. self._sd.nFluids .. " elements.")
   assert(#tbl.mass == self._sd.nFluids, "Updater.FiveMomentSrc: Mass table must have " .. self._sd.nFluids .. " elements.")

   self._sd.epsilon0 = assert(tbl.epsilon0, "Updater.FiveMomentSrc: Must specify 'epsilon0'")
   self._sd.chi_e = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self._sd.chi_m = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0

   self._sd.gravity = tbl.gravity and tbl.gravity or 0.0
   self._sd.gravityDir = tbl.dir and tbl.dir or 1 -- by default gravity acts in X direction
   
   self._sd.hasStatic = tbl.hasStaticField ~= nil and tbl.hasStaticField or false
   self._sd.hasPressure = tbl.hasPressure ~= nil and tbl.hasPressure or true

   self._fd = ffi.new("FluidData_t[?]", self._sd.nFluids)
   -- store charge and mass for each fluid
   for n = 1, self._sd.nFluids do
      self._fd[n-1].charge = tbl.charge[n]
      self._fd[n-1].mass = tbl.mass[n]
   end

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
end

-- advance method
function FiveMomentSrc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local nFluids = self._sd.nFluids
   
   -- check if correct number of inputs were provided
   assert(#outFld == nFluids+1,
	  "Must have exactly " .. nFluids+1 .. " output fields. Provided " .. #inFld .. " instead.")
   local emFld = outFld[#outFld] -- EM field

   -- make list of indexer functions
   local fIdxr = {}
   for i = 1, nFluids do
      fIdxr[i] = outFld[i]:genIndexer()
   end
   local emIdxr = emFld:genIndexer()

   -- allocate stuff to pass to C
   local fDp = ffi.new("double*[?]", nFluids)
   local emDp = ffi.new("double*")

   local localRange = emFld:localRange()   
   -- loop over local range, updating source in each cell
   for idx in localRange:colMajorIter() do
      -- set pointers to fluids and field
      for i = 1, nFluids do
	 fDp[i-1] = outFld[i]:getDataPtrAt(fIdxr[i](idx))
      end
      emDp = emFld:getDataPtrAt(emIdxr(idx))

      -- update sources
      self._updateSrc(self, dt, fDp, emDp)
   end

   return true, GKYL_MAX_DOUBLE
end

return FiveMomentSrc
