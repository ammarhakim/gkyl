--------------------------------------------------------------------------------
-- Updater to apply inter-species friction for the multifluid five-moment
-- (Euler) model.
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
    double mass;
  } FiveMomentFrictionFluidData_t;

  typedef struct {
    int nFluids;
    double gasGamma;
    bool hasPressure;
    double nuBase[5*4/2];
  } FiveMomentFrictionSrcData_t;

  void gkylFiveMomentFrictionSrcForwardEuler(
    const FiveMomentFrictionSrcData_t *sd,
    const FiveMomentFrictionFluidData_t *fd,
    const double dt,
    double **fPtrs);

  void gkylFiveMomentFrictionSrcTimeCentered(
    const FiveMomentFrictionSrcData_t *sd,
    const FiveMomentFrictionFluidData_t *fd,
    const double dt,
    double **fPtrs);

  void gkylFiveMomentFrictionSrcExact(
    const FiveMomentFrictionSrcData_t *sd,
    const FiveMomentFrictionFluidData_t *fd,
    const double dt,
    double **fPtrs);
]]

local function forwardEuler(self, dt, fPtrs)
   ffi.C.gkylFiveMomentFrictionSrcForwardEuler(self._sd, self._fd, dt, fPtrs)
end

local function timeCentered(self, dt, fPtrs)
   ffi.C.gkylFiveMomentFrictionSrcTimeCentered(self._sd, self._fd, dt, fPtrs)
end

local function exact(self, dt, fPtrs)
   ffi.C.gkylFiveMomentFrictionSrcExact(self._sd, self._fd, dt, fPtrs)
end

local FiveMomentFrictionSrc = Proto(UpdaterBase)

function FiveMomentFrictionSrc:init(tbl)
   FiveMomentFrictionSrc.super.init(self, tbl)

   local pfx = "Updater.FiveMomentFrictionSrc: "

   self._onGrid = assert(tbl.onGrid, pfx.."Must set grid object bt 'onGrid'.")

   local nFluids = assert(tbl.numFluids, pfx.."Must specify 'numFluids'.")
   self._sd = ffi.new(ffi.typeof("FiveMomentFrictionSrcData_t"))
   self._sd.nFluids = nFluids
   self._sd.gasGamma = assert(tbl.gasGamma, pfx.."Must specify 'gasGamma'.")
   self._sd.hasPressure = xsys.pickBool(tbl.hasPressure, true)
  
   local nuBase = assert(tbl.nu, pfx.."Must specify 'nu' table.")
   assert(#nuBase==nFluids*(nFluids-1)/2, pfx.."'nu' entry # is incorrect.")
   -- FIXME: Presently we are statically allocating space for _sd.nuBase (see
   -- definition of FiveMomentFrictionSrcData_t); Dynamically allocating nuBase
   -- using ffi.new within the dynamically allocated _sd seems to possibly cause
   -- memory corruption.
   if false then
      self._sd.nuBase = ffi.new("double[?]", #nuBase)
   else
      assert(nFluids<=5, pfx.."Presently up to 5 species are supported.")
   end
   for i=1,#nuBase do
      self._sd.nuBase[i-1] = nuBase[i]
   end

   assert(#tbl.mass == nFluids,
          pfx.."mass table must have "..nFluids.." entries.")
   self._fd = ffi.new("FiveMomentFrictionFluidData_t[?]", nFluids)
   for s = 1, nFluids do
      self._fd[s-1].mass = tbl.mass[s]
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
   elseif scheme == "time-centered" then
      self._updateSrc = timeCentered
   elseif scheme == "exact" then
      assert(nFluids==2, pfx.."For the 'exact' scheme, nFluids must be 2.")
      self._updateSrc = exact
   else
      assert(false, string.format(pfx.."Scheme %s is not supported", scheme))
   end
   self.scheme = scheme
end

function FiveMomentFrictionSrc:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local ndim = grid:ndim()
   local nFluids = #outFld
   assert(nFluids == self._sd.nFluids)

   local fIdxrs = {}
   for i = 1, nFluids do
      fIdxrs[i] = outFld[i]:genIndexer()
   end
   local fPtrs = ffi.new("double*[?]", nFluids)

   for idx in outFld[1]:localRangeIter() do
      grid:setIndex(idx)

      for i = 1, nFluids do
         fPtrs[i-1] = outFld[i]:getDataPtrAt(fIdxrs[i](idx))
      end

      self._updateSrc(self, dt, fPtrs)
   end

   return true, GKYL_MAX_DOUBLE
end

return FiveMomentFrictionSrc
