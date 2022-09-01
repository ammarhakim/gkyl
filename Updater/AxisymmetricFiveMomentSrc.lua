-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update geometric sources for the five-moment (i.e., Euler) model
-- to effectively achieve axisymmetry.
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
    bool evolve;
  } AxisymmetricFluidData_t;

  typedef struct {
    int nFluids;
    double gasGamma;
    bool hasPressure;
  } AxisymmetricFiveMomentSrcData_t;

  void gkylAxisymmetricFiveMomentSrcForwardEuler(
    const AxisymmetricFiveMomentSrcData_t *sd,
    const AxisymmetricFluidData_t *fd,
    const double dt,
    const double *xc,
    double **fPtrs);

  void
  gkylAxisymmetricFiveMomentSrcRk3(const AxisymmetricFiveMomentSrcData_t *sd,
                                   const AxisymmetricFluidData_t *fd,
                                   const double dt,
                                   const double *xc,
                                   double **fPtrs);

  void
  gkylAxisymmetricFiveMomentSrcSemiExact(
    const AxisymmetricFiveMomentSrcData_t *sd,
    const AxisymmetricFluidData_t *fd,
    const double dt,
    const double *xc,
    double **fPtrs);
]]

local function forwardEuler(self, dt, xc, fPtrs)
   ffi.C.gkylAxisymmetricFiveMomentSrcForwardEuler(self._sd, self._fd, dt, xc, fPtrs)
end

local function rk3(self, dt, xc, fPtrs)
   ffi.C.gkylAxisymmetricFiveMomentSrcRk3(self._sd, self._fd, dt, xc, fPtrs)
end

local function semiExact(self, dt, xc, fPtrs)
   ffi.C.gkylAxisymmetricFiveMomentSrcSemiExact(self._sd, self._fd, dt, xc, fPtrs)
end

local AxisymmetricFiveMomentSrc = Proto(UpdaterBase)

function AxisymmetricFiveMomentSrc:init(tbl)
   AxisymmetricFiveMomentSrc.super.init(self, tbl)

   self._onGrid = assert(tbl.onGrid, "Updater.AxisymmetricFiveMomentSrc: Must provide grid object using 'onGrid'")

   self._sd = ffi.new(ffi.typeof("AxisymmetricFiveMomentSrcData_t"))
   self._sd.nFluids = assert(tbl.numFluids, "Updater.AxisymmetricFiveMomentSrc: Must specify number of fluid using 'numFluids'")
   self._sd.gasGamma = assert(tbl.gasGamma, "Updater.AxisymmetricFiveMomentSrc: Must specify 'gasGamma'")
   self._sd.hasPressure = tbl.hasPressure ~= nil and tbl.hasPressure or true

   self._fd = ffi.new("AxisymmetricFluidData_t[?]", self._sd.nFluids)
   for s = 1, self._sd.nFluids do
      if (tbl.evolve ~= nil) then
        self._fd[s-1].evolve = tbl.evolve[s]
      else
        self._fd[s-1].evolve = true
      end
   end

   local scheme = tbl.scheme and tbl.scheme or "forwardEuler"
   self._updateSrc = nil
   if scheme == "forwardEuler" then
      self._updateSrc = forwardEuler
   elseif scheme == "rk3" then
      self._updateSrc = rk3
   elseif scheme == "semi-exact" then
      self._updateSrc = semiExact
   else
      assert(false, string.format("Updater.AxisymmetricFiveMomentSrc: scheme %s not supported", scheme))
   end
   self.scheme = scheme
end

function AxisymmetricFiveMomentSrc:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local ndim = grid:ndim()
   local xc = Lin.Vec(ndim)
   local nFluids = #outFld
   assert(nFluids == self._sd.nFluids)

   local fIdxrs = {}
   for i = 1, nFluids do
      fIdxrs[i] = outFld[i]:genIndexer()
   end
   local fPtrs = ffi.new("double*[?]", nFluids)

   for idx in outFld[1]:localRangeIter() do
      grid:setIndex(idx)
      grid:cellCenter(xc)

      for i = 1, nFluids do
         fPtrs[i-1] = outFld[i]:getDataPtrAt(fIdxrs[i](idx))
      end

      self._updateSrc(self, dt, xc:data(), fPtrs)
   end

   return self:reduceStatusDt(true, GKYL_MAX_DOUBLE)
end

return AxisymmetricFiveMomentSrc
