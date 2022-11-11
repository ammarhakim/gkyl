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
} FiveMomentUserSrcFluidData_t;

typedef struct {
 int nFluids;
 double gasGamma;
} FiveMomentUserSrcData_t;

void gkylFiveMomentUserSrcForwardEuler(
  const FiveMomentUserSrcData_t *sd,
  const FiveMomentUserSrcFluidData_t *fd,
  const double dt,
  double *qPtrs[],
  const double *sourcePtrs[]);
]]

local function forwardEuler(self, dt, fPtrs, sourcePtrs)
   ffi.C.gkylFiveMomentUserSrcForwardEuler(
      self._sd, self._fd, dt, fPtrs, sourcePtrs)
end

local FiveMomentUserSrc = Proto(UpdaterBase)

function FiveMomentUserSrc:init(tbl)
   FiveMomentUserSrc.super.init(self, tbl)

   local pfx = "Updater.FiveMomentUserSrc: "

   self._onGrid = assert(tbl.onGrid, pfx.."Must set grid object bt 'onGrid'.")

   local nFluids = assert(tbl.numFluids, pfx.."Must specify 'numFluids'.")
   self._sd = ffi.new(ffi.typeof("FiveMomentUserSrcData_t"))
   self._sd.nFluids = nFluids
   self._sd.gasGamma = assert(tbl.gasGamma, pfx.."Must specify 'gasGamma'.")

   assert(#tbl.mass == nFluids,
          pfx.."mass table must have "..nFluids.." entries.")
   self._fd = ffi.new("FiveMomentUserSrcFluidData_t[?]", nFluids)
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
   else
      assert(false, string.format(pfx.."Scheme %s is not supported", scheme))
   end
   self.scheme = scheme
end

function FiveMomentUserSrc:_advance(tCurr, srcFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local ndim = grid:ndim()
   local nFluids = #outFld
   assert(nFluids == self._sd.nFluids)

   local fIdxrs, sIdxrs = {}, {}
   for i = 1, nFluids do
      fIdxrs[i] = outFld[i]:genIndexer()
      sIdxrs[i] = srcFld[i]:genIndexer()
   end
   local fPtrs = ffi.new("double*[?]", nFluids)
   local sPtrs = ffi.new("const double*[?]", nFluids)

   for idx in outFld[1]:localRangeIter() do
      grid:setIndex(idx)

      for i = 1, nFluids do
         fPtrs[i-1] = outFld[i]:getDataPtrAt(fIdxrs[i](idx))
         sPtrs[i-1] = srcFld[i]:getDataPtrAt(sIdxrs[i](idx))
      end

      self._updateSrc(self, dt, fPtrs, sPtrs)
   end

   return true, GKYL_MAX_DOUBLE
end

return FiveMomentUserSrc
