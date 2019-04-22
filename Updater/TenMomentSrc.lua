-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update ten-moment source terms. This updater allows
-- both explicit and implicit updates.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local COL_PIV_HOUSEHOLDER_QR = 0;
local PARTIAL_PIV_LU = 1;

-- system libraries
local ffi = require "ffi"

-- Define C types for storing private data for use in updater
ffi.cdef [[
typedef struct {
  double charge, mass; /* Charge and mass */
  bool evolve;

  double qbym;
} FluidData_t;

typedef struct {
  int8_t nFluids; /* Number of fluids */
  double epsilon0; /* Permittivity of free space */
  double chi_e, chi_m; /* Propagation speed factor for electric field error potential */
  int8_t gravityDir; /* Direction of gravity force */
  double gravity; /* Gravitational acceleration */
  bool hasStatic, hasPressure; /* Flag to indicate if there is: static EB field, pressure */
  bool hasSigma; /* Flag to indicate if there is: sigma */
  int8_t linSolType; /* Flag to indicate linear solver type for implicit method */
} MomentSrcData_t;

  void gkylTenMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma);
  void gkylTenMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma);
  void gkylTenMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma);
  void gkylTenMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma);
  void gkylTenMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma);
]]

-- Explicit, SSP RK3 scheme
local function updateSrcRk3(self, dt, fPtr, emPtr)
   ffi.C.gkylTenMomentSrcRk3(self._sd, self._fd, dt, fPtr, emPtr)
end

-- Use an explicit scheme to update momentum and electric field: this
-- is an extension of the standard Boris push algorithm, in which the
-- half time-step electric field update is replaced by an implicit
-- step in which both the velocity and electric are updated.
local function updateSrcModBoris(self, dt, fPtr, emPtr)
   print("updateSrcModBoris")
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCentered(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
   ffi.C.gkylTenMomentSrcTimeCentered(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCenteredDirect2(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
   ffi.C.gkylTenMomentSrcTimeCenteredDirect2(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCenteredDirect(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
   ffi.C.gkylTenMomentSrcTimeCenteredDirect(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
end

-- Use an exact scheme to update momentum and electric field
local function updateSrcExact(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
   ffi.C.gkylTenMomentSrcExact(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr)
end

-- Ten-moment source updater object
local TenMomentSrc = Proto(UpdaterBase)

function TenMomentSrc:init(tbl)
   TenMomentSrc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.TenMomentSrc: Must provide grid object using 'onGrid'")

   self._sd = ffi.new(ffi.typeof("MomentSrcData_t"))   
   -- Read in solver parameters
   self._sd.nFluids = assert(tbl.numFluids, "Updater.TenMomentSrc: Must specify number of fluid using 'numFluids'")

   assert(#tbl.charge == self._sd.nFluids, "Updater.TenMomentSrc: Charge table must have " .. self._sd.nFluids .. " elements.")
   assert(#tbl.mass == self._sd.nFluids, "Updater.TenMomentSrc: Mass table must have " .. self._sd.nFluids .. " elements.")

   self._sd.epsilon0 = assert(tbl.epsilon0, "Updater.TenMomentSrc: Must specify 'epsilon0'")
   self._sd.chi_e = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self._sd.chi_m = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 1.0

   self._sd.gravity = tbl.gravity and tbl.gravity or 0.0
   self._sd.gravityDir = tbl.dir and tbl.dir or 1 -- by default gravity acts in X direction
   
   self._sd.hasStatic = tbl.hasStaticField ~= nil and tbl.hasStaticField or false
   self._sd.hasPressure = tbl.hasPressure ~= nil and tbl.hasPressure or true
   self._sd.hasSigma = tbl.hasSigmaField ~= nil and tbl.hasSigmaField or false
   
   self._sigmaFld = tbl.sigmaField
   
   local linSolType = tbl.linSolType and tbl.linSolType or "partialPivLu"
   if linSolType == "partialPivLu" then
      self._sd.linSolType = PARTIAL_PIV_LU
   elseif linSolType == "colPivHouseholderQr" then
      self._sd.linSolType = COL_PIV_HOUSEHOLDER_QR
   else
     assert(false, string.format("linSolType %s not supported", linSolType))
   end

   self._fd = ffi.new("FluidData_t[?]", self._sd.nFluids)
   -- store charge and mass for each fluid
   for n = 1, self._sd.nFluids do
      self._fd[n-1].charge = tbl.charge[n]
      self._fd[n-1].mass = tbl.mass[n]
      self._fd[n-1].qbym = self._fd[n-1].charge / self._fd[n-1].mass
      -- self._fd[n-1].evolve = tbl.evolve ~= nil and tbl.evolve[n] or true
      if (tbl.evolve ~= nil) then
        self._fd[n-1].evolve = tbl.evolve[n]
      else
        self._fd[n-1].evolve = true
      end
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
   elseif scheme == "time-centered-direct2" or scheme == "direct2" then
      self._updateSrc = updateSrcTimeCenteredDirect2
   elseif scheme == "time-centered-direct" or scheme == "direct" then
      self._updateSrc = updateSrcTimeCenteredDirect
   elseif scheme == "exact" then
      self._updateSrc = updateSrcExact
   end
end

-- advance method
function TenMomentSrc:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local nFluids = self._sd.nFluids
   
   local ndim = grid:ndim()
   local xc = Lin.Vec(ndim) -- cell center
   
   -- check if correct number of inputs were provided
   assert(#outFld == nFluids+1,
	  "Must have exactly " .. nFluids+1 .. " output fields. Provided " .. #outFld .. " instead.")
   local emFld = outFld[#outFld] -- EM field

   local staticEmFld
   if (self._sd.hasStatic) then
      staticEmFld = inFld[1] -- static EM field
   end

   -- make list of indexer functions
   local fIdxr = {}
   for i = 1, nFluids do
      fIdxr[i] = outFld[i]:genIndexer()
   end
   local emIdxr = emFld:genIndexer()
   local staticEmIdxr
   if (self._sd.hasStatic) then
      staticEmIdxr = staticEmFld:genIndexer()
   end

   local sigmaIdxr
   if (self._sd.hasSigma) then
      sigmaIdxr = self._sigmaFld:genIndexer()
   end

   -- allocate stuff to pass to C
   local fDp = ffi.new("double*[?]", nFluids)
   local emDp = ffi.new("double*")
   local staticEmDp = ffi.new("double*")
   local sigmaDp = ffi.new("double*")

   local localRange = emFld:localRange()   
   -- loop over local range, updating source in each cell
   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)
      grid:cellCenter(xc)

      -- set pointers to fluids and field
      for i = 1, nFluids do
	 fDp[i-1] = outFld[i]:getDataPtrAt(fIdxr[i](idx))
      end
      emDp = emFld:getDataPtrAt(emIdxr(idx))
      if (self._sd.hasStatic) then
         staticEmDp = emFld:getDataPtrAt(staticEmIdxr(idx))
      end
      if (self._sd.hasSigma) then
         sigmaDp = self._sigmaFld:getDataPtrAt(sigmaIdxr(idx))
      end

      -- update sources
      self._updateSrc(self, dt, fDp, emDp, staticEmDp, sigmaDp)
   end

   return true, GKYL_MAX_DOUBLE
end

return TenMomentSrc
