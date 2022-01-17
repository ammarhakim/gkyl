-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute RHS or forward Euler update for hyperbolic
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local DataStruct = require "DataStruct"
local Grid = require "Grid.RectCart"
local CartField = require "DataStruct.CartField"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local cuda = nil
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
   cuAlloc = require "Cuda.Alloc"
end

ffi.cdef [[ 
  typedef struct GkylEquation_t GkylEquation_t ;
  typedef struct {
      int updateDirs[6];
      bool zeroFluxFlags[6];
      int32_t numUpdateDirs;
      bool updateVolumeTerm;
      double dt;
      GkylEquation_t *equation;
      GkylCartField_t *cflRateByCell;
      GkylCartField_t *maxsByCell;
      double *maxs;
  } GkylHyperDisCont_t; 

  void advanceOnDevice(const int numBlocks, const int numThreads, const int numComponents, const GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut);
  void advanceOnDevice_shared(int numBlocks, int numThreads, int numComponents, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut);
  void setDtAndCflRate(GkylHyperDisCont_t *hyper, double dt, GkylCartField_t *cflRate);
]]

-- Hyperbolic DG solver updater object
local HyperDisCont = Proto(UpdaterBase)

function HyperDisCont:init(tbl)
   HyperDisCont.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.HyperDisCont: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.HyperDisCont: Must specify basis functions to use using 'basis'")

   -- By default, clear output field before incrementing with vol/surf updates.
   self._clearOut = xsys.pickBool(tbl.clearOut, true)

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._onGrid:ndim()

   -- Equation to solve.
   self._equation = assert(tbl.equation, "Updater.HyperDisCont: Must provide equation object using 'equation'")

   -- By default, update all directions.
   self._updateDirs = {}
   for d = 1, self._ndim do
      self._updateDirs[d] = d
   end
   -- Read in which directions we are to update.
   if tbl.updateDirections then
      self._updateDirs = tbl.updateDirections
   end

   -- Set zero flux direction flags.
   self._zeroFluxFlags = {}
   for d = 1, self._ndim do self._zeroFluxFlags[d] = false end
   if tbl.zeroFluxDirections then
      for i, d in ipairs(tbl.zeroFluxDirections) do self._zeroFluxFlags[d] = true end
   end

   -- Flag to turn on/off volume term.
   self._updateVolumeTerm = xsys.pickBool(tbl.updateVolumeTerm, true)

   -- Flag to indicate the use of local or global upwind fluxes. Local
   -- does not need to do a reduction of the maximum speed (saves an Mpi.Allreduce).
   self._globalUpwind = xsys.pickBool(tbl.globalUpwind, true)

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.HyperDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- Maximum characteristic velocities for use in penalty based fluxes.
   self._maxs, self._maxsOld, self._maxsLocal = Lin.Vec(self._ndim), Lin.Vec(self._ndim), Lin.Vec(self._ndim)
   for d = 1, self._ndim do
      -- Very first step the penalty term will be zero. However, in
      -- subsequent steps the maximum speed from the previous step
      -- will be used.
      self._maxs[d] = 0.0
      self._maxsOld[d] = 0.0
   end
   self._noPenaltyFlux = xsys.pickBool(tbl.noPenaltyFlux, false)

   -- Mask field, M=1 for cells where updates should be computed,
   -- and M=0 in cells were updates are not to be calculated.
   self._maskFld = tbl.maskField
   self._maskM, self._maskP = self._maskFld:get(1), self._maskFld:get(1)
   self._indxrP0 = self._maskFld:genIndexer()

   self._isFirst = true
   self._auxFields = {} -- Auxilliary fields passed to eqn object.
   self._perpRangeDecomp = {} -- Perp ranges in each direction.

   -- to store grid info (pre-allocated to avoid allocations in time loop).
   self._dxp,  self._dxm  = Lin.Vec(self._ndim),    Lin.Vec(self._ndim)    -- cell shape on right/left
   self._xcp,  self._xcm  = Lin.Vec(self._ndim),    Lin.Vec(self._ndim)    -- cell center on right/left
   self._idxp, self._idxm = Lin.IntVec(self._ndim), Lin.IntVec(self._ndim) -- index on right/left


   return self
end

function HyperDisCont:initDevice(tbl)
   self.maxsByCell = DataStruct.Field {
      onGrid = self._onGrid,
      numComponents = self._ndim,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   self.maxs = cuAlloc.Double(self._ndim)
   local hyper = ffi.new("GkylHyperDisCont_t")
   hyper.updateDirs = ffi.new("int[6]", self._updateDirs)
   hyper.zeroFluxFlags = ffi.new("bool[6]", self._zeroFluxFlags)
   hyper.numUpdateDirs = #self._updateDirs
   hyper.updateVolumeTerm = self._updateVolumeTerm
   hyper.equation = self._equation._onDevice
   hyper.maxsByCell = self.maxsByCell._onDevice
   hyper.maxs = self.maxs:data()
   self._onHost = hyper
   local sz = sizeof("GkylHyperDisCont_t")
   self._onDevice, err = cuda.Malloc(sz)
   cuda.Memcpy(self._onDevice, hyper, sz, cuda.MemcpyHostToDevice)

   self.numThreads = tbl.numThreads or GKYL_DEFAULT_NUM_THREADS
   self._useSharedDevice = xsys.pickBool(tbl.useSharedDevice, false)

   return self
end

local isInside = function (inOutPtr) return inOutPtr:data()[0] > 0. end

-- advance method
function HyperDisCont:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local cflRateByCell = self._cflRateByCell

   local qIn = assert(inFld[1], "HyperDisCont.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "HyperDisCont.advance: Must specify an output field")

   -- pass aux fields to equation object
   for i = 1, #inFld-1 do self._auxFields[i] = inFld[i+1] end
   self._equation:setAuxFields(self._auxFields)

   local ndim = grid:ndim()

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local localRange = qRhsOut:localRange()
   local globalRange = qRhsOut:globalRange()
   local indxr = qIn:genIndexer() -- indexer functions into fields
   local indxrP0 = self._indxrP0

   -- to store grid info
   local dxp,  dxm  = self._dxp,  self._dxm  -- cell shape on right/left
   local xcp,  xcm  = self._xcp,  self._xcm  -- cell center on right/left
   local idxp, idxm = self._idxp, self._idxm -- index on right/left

   -- pointers for (re)use in update
   local qInM, qInP = qIn:get(1), qIn:get(1)
   local qRhsOutM, qRhsOutP = qRhsOut:get(1), qRhsOut:get(1)
   local cflRateByCellP = cflRateByCell:get(1)
   local cflRateByCellM = cflRateByCell:get(1)

   -- This flag is needed as the volume integral already contains
   -- contributions from all directions. Hence, we must only
   -- accumulate the volume contribution once, skipping it for other
   -- directions.
   local firstDir = true

   -- Use maximum characteristic speeds from previous step as penalty.
   for d = 1, ndim do
      if not self._noPenaltyFlux then self._maxsOld[d] = self._maxs[d] end
      self._maxsLocal[d] = 0.0 -- Reset to get new values in this step.
   end

   local tId = grid:subGridSharedId() -- Local thread ID.

   -- Clear output field before computing vol/surf increments.
   if self._clearOut then qRhsOut:clear(0.0) end
   -- Accumulate contributions from volume and surface integrals.
   local cflRate
   local isInsideM, isInsideP  -- Flags indicating if cells are inside the domain to be updated.
   -- Iterate through updateDirs backwards so that a zero flux dir is first in kinetics.
   for i = #self._updateDirs, 1, -1 do 
      local dir = self._updateDirs[i]
      -- Lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell).
      local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
      local dirLoSurfIdx, dirUpSurfIdx = dirLoIdx, dirUpIdx
      
      -- Compute loop bounds for zero flux direction.
      if self._zeroFluxFlags[dir] then
         local dirGlobalLoIdx, dirGlobalUpIdx = globalRange:lower(dir), globalRange:upper(dir)+1
         if dirLoIdx == dirGlobalLoIdx then dirLoSurfIdx = dirLoIdx+1 end
         if dirUpIdx == dirGlobalUpIdx then dirUpSurfIdx = dirUpIdx-1 end
      end

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
      for idx in perpRangeDecomp:rowMajorIter(tId) do
	 idx:copyInto(idxp); idx:copyInto(idxm)

         for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    idxm[dir], idxp[dir] = i-1, i -- cell left/right of edge 'i'

	    grid:setIndex(idxm)
            grid:getDx(dxm)
	    grid:cellCenter(xcm)

	    grid:setIndex(idxp)
            grid:getDx(dxp)
	    grid:cellCenter(xcp)

	    qIn:fill(indxr(idxm), qInM)
	    qIn:fill(indxr(idxp), qInP)

	    qRhsOut:fill(indxr(idxm), qRhsOutM)
	    qRhsOut:fill(indxr(idxp), qRhsOutP)
            cflRateByCell:fill(indxr(idxm), cflRateByCellM)
            cflRateByCell:fill(indxr(idxp), cflRateByCellP)

	    self._maskFld:fill(indxrP0(idxm), self._maskM)
	    self._maskFld:fill(indxrP0(idxp), self._maskP)
            isInsideM, isInsideP = isInside(self._maskM), isInside(self._maskP)

	    if firstDir and i<=dirUpIdx-1 and self._updateVolumeTerm and isInsideP then
	       cflRate = self._equation:volTerm(xcp, dxp, idxp, qInP, qRhsOutP)
               cflRateByCellP:data()[0] = cflRateByCellP:data()[0] + cflRate
	    end

	    if i >= dirLoSurfIdx and i <= dirUpSurfIdx and (isInsideM and isInsideP) then
               local cflp = cflRateByCellP:data()[0]*dt/.9 -- .9 here is conservative, but we are using dt from the prev step
               local cflm = cflRateByCellM:data()[0]*dt/.9 -- .9 here is conservative, but we are using dt from the prev step
               if cflp == 0.0 then cflp = math.min(1.5*cflm, cfl) end
               if cflm == 0.0 then cflm = math.min(1.5*cflp, cfl) end
	       local maxs = self._equation:surfTerm(
	          dir, cflm, cflp, xcm, xcp, dxm, dxp, self._maxsOld[dir], idxm, idxp, qInM, qInP, qRhsOutM, qRhsOutP)
	       self._maxsLocal[dir] = math.max(self._maxsLocal[dir], maxs)
            else
	       if self._zeroFluxFlags[dir] and (isInsideP or isInsideM) then
	          -- we need to give equations a chance to apply partial surface
	          -- updates even when the zeroFlux BCs have been applied
                  local edge
                  if self._maskFld then
                     edge = (isInsideP and not isInsideM) and -1 or 1
                  else
                     edge = i<dirLoSurfIdx and -1 or 1
                  end
	          self._equation:boundarySurfTerm(
	             dir, xcm, xcp, dxm, dxp, self._maxsOld[dir], idxm, idxp, edge, qInM, qInP, qRhsOutM, qRhsOutP)
               end
	    end
	 end
      end
      if firstDir then cflRateByCell:sync() end
      firstDir = false
   end

   -- Determine largest amax across processors.
   if self._globalUpwind then
      Mpi.Allreduce(
         self._maxsLocal:data(), self._maxs:data(), ndim, Mpi.DOUBLE, Mpi.MAX, self:getComm())
   end

   self._isFirst = false
end

function HyperDisCont:_advanceOnDevice(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "HyperDisCont.advanceOnDevice: Must specify an input field")
   local qRhsOut = assert(outFld[1], "HyperDisCont.advanceOnDevice: Must specify an output field")

   for i = 1, #inFld-1 do
      self._auxFields[i] = inFld[i+1]
   end

   self._equation:setAuxFieldsOnDevice(self._auxFields)

   local numCellsLocal = qRhsOut:localRange():volume()
   local numThreads = math.min(self.numThreads, numCellsLocal)
   local numBlocks  = math.ceil(numCellsLocal/numThreads)

   if self._clearOut then
      cuda.Memset(qRhsOut:deviceDataPointer(), 0.0, sizeof('double')*qRhsOut:size())
   end

   if self._useSharedDevice then
      ffiC.advanceOnDevice_shared(numBlocks, numThreads, qIn:numComponents(), self._onDevice, qIn._onDevice, qRhsOut._onDevice)
   else
      ffiC.advanceOnDevice(numBlocks, numThreads, qIn:numComponents(), self._onDevice, qIn._onDevice, qRhsOut._onDevice)
   end

   self.maxsByCell:deviceReduce('max', self.maxs)  
end

-- set up pointers to dt and cflRateByCell
function HyperDisCont:setDtAndCflRate(dt, cflRateByCell)
   HyperDisCont.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

return HyperDisCont
