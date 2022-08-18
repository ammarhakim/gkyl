-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update geometric sources for the perfectly-hyperbolic Maxwell's
-- equations to effectively achieve axisymmetry.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype =
  xsys.from(ffi, "new, copy, fill, sizeof, typeof, metatype")

ffi.cdef [[
  typedef struct {
    double epsilon0;
    double mu0;
    double chi_e;
    double chi_m;
  } AxisymmetricPhMaxwellSrcData_t;

  void gkylAxisymmetricPhMaxwellSrcForwardEuler(
    const AxisymmetricPhMaxwellSrcData_t *sd,
    const double dt,
    const double *xc,
    double *fPtr);

  void
  gkylAxisymmetricPhMaxwellSrcRk3(const AxisymmetricPhMaxwellSrcData_t *sd,
                                  const double dt,
                                  const double *xc,
                                  double *fPtr);

]]

local function forwardEuler(self, dt, xc, fPtr)
   ffi.C.gkylAxisymmetricPhMaxwellSrcForwardEuler(self._sd, dt, xc, fPtr)
end

local function rk3(self, dt, xc, fPtr)
   ffi.C.gkylAxisymmetricPhMaxwellSrcRk3(self._sd, dt, xc, fPtr)
end

local AxisymmetricPhMaxwellSrc = Proto(UpdaterBase)

function AxisymmetricPhMaxwellSrc:init(tbl)
   AxisymmetricPhMaxwellSrc.super.init(self, tbl)

   self._onGrid = assert(tbl.onGrid, "Updater.AxisymmetricPhMaxwellSrc: Must provide grid object using 'onGrid'")

   self._sd = ffi.new(ffi.typeof("AxisymmetricPhMaxwellSrcData_t"))   
   self._sd.epsilon0 = assert(tbl.epsilon0, "Updater.AxisymmetricPhMaxwellSrc: Must specify 'epsilon0'")
   self._sd.mu0 = assert(tbl.mu0, "Updater.AxisymmetricPhMaxwellSrc: Must specify 'mu0'")
   self._sd.chi_e = tbl.chi_e and tbl.chi_e or 0
   self._sd.chi_m = tbl.chi_m and tbl.chi_m or 1

   local scheme = tbl.scheme and tbl.scheme or "forwardEuler"
   self._updateSrc = nil
   if scheme == "forwardEuler" then
      self._updateSrc = forwardEuler
   elseif scheme == "rk3" then
      self._updateSrc = rk3
   else
      assert(false, string.format("Updater.AxisymmetricPhMaxwellSrc: scheme %s not supported", scheme))
   end
   self.scheme = scheme
end

function AxisymmetricPhMaxwellSrc:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local ndim = grid:ndim()
   local xc = Lin.Vec(ndim)
   local emf = outFld[1]

   local fIndexer = emf:genIndexer()
   -- local fPtr = ffi.new("double*")

   for idx in emf:localRangeIter() do
      grid:setIndex(idx)
      grid:cellCenter(xc)
      local fPtr = emf:getDataPtrAt(fIndexer(idx))

      self._updateSrc(self, dt, xc:data(), fPtr)
   end

   return self:reduceStatusDt(true, GKYL_MAX_DOUBLE)
end

return AxisymmetricPhMaxwellSrc
