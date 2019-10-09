-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions at stair-stepped boundaries.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries
local Alloc        = require "Lib.Alloc"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"

-- System libraries.
local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
"new, copy, fill, sizeof, typeof, metatype")

-- Helper object for indexing 1D slice data. The slice spans from
-- [lower, upper] (inclusive) and has `stride` pieces of data stored
-- at each location.
local createSliceData = function (dtype)
   local slice_mt = {
      __new = function (self, lower, upper, stride)
         local n = upper-lower+1
         local v = new(self)
         v._data = Alloc.malloc(typeof(dtype)*n*stride)
         v._sz, v._stride, v._lower = n*stride, stride, lower
         return v
      end,
      __index = function (self, k)
         return self._data+(k-self._lower)*self._stride
      end,
      __gc = function (self)
         Alloc.free(self._data)
      end,
   }
   return metatype(typeof(string.format("struct {int32_t _sz, _stride, _lower; %s *_data; }", dtype)), slice_mt)
end
local SliceDataBool = createSliceData("bool")

-- Boundary condition updater
local StairSteppedBc = Proto(UpdaterBase)

local isOutside = function (inOutPtr)
   return inOutPtr[1] < 0
end

function StairSteppedBc:init(tbl)
   StairSteppedBc.super.init(self, tbl) -- setup base object

   self._isFirst = true -- will be reset first time _advance() is called

   self._grid = assert(tbl.onGrid, "Updater.StairSteppedBc: Must specify grid to use with 'onGrid''")
   self._dir  = assert(tbl.dir, "Updater.StairSteppedBc: Must specify direction to apply BCs with 'dir'")

   self._bcList = assert(
   tbl.boundaryConditions, "Updater.StairSteppedBc: Must specify boundary conditions to apply with 'boundaryConditions'")

   self._localRange = tbl.onGrid:localRange()
   local localRange = self._localRange

   self._inOut = assert(tbl.inOut, "Updater.StairSteppedBc: Must specify mask field with 'inOut'")

   local l, u = localRange:lower(self._dir)-2, localRange:upper(self._dir)+2
   self._isGhostL = SliceDataBool(l, u, 1);
   self._isGhostR = SliceDataBool(l, u, 1);
end

function StairSteppedBc:_advance(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "StairSteppedBc.advance: Must-specify an output field")
   
   local dir = self._dir
   local localRange = grid:localRange()

   if self._isFirst then
      self._perpRangeDecomp = LinearDecomp.LinearDecompRange {
         range      = localRange:shorten(dir), -- range orthogonal to 'dir'
         numSplit   = grid:numSharedProcs(),
         threadComm = self:getSharedComm()
      }
   end

   -- get pointers to (re)use inside inner loop [G: Ghost, S: Skin]
   local qG, qS  = qOut:get(1), qOut:get(1)
   local idxL    = Lin.IntVec(grid:ndim())
   local idxR    = Lin.IntVec(grid:ndim())
   local indexer = qOut:genIndexer()

   local inOut          = self._inOut
   local inOutL, inOutR = inOut:get(1), inOut:get(1)
   local inOutIdxr      = inOut:genIndexer()
   local isGhostL       = self._isGhostL
   local isGhostR       = self._isGhostR

   local tId = self._grid:subGridSharedId() -- local thread ID

   -- outer loop is over directions orthogonal to 'dir' and inner
   -- loop is over 1D slice in `dir`.
   for idx in self._perpRangeDecomp:colMajorIter(tId) do
      idx:copyInto(idxL)
      idx:copyInto(idxR)

      -- loop over edges
      for i = localRange:lower(dir)-1, localRange:upper(dir)+1 do
         idxL[dir] = i
         idxR[dir] = i + 1
         inOut:fill(inOutIdxr(idxL), inOutL)
         inOut:fill(inOutIdxr(idxR), inOutR)
         local isOutsideL = isOutside(inOutL)
         local isOutsideR = isOutside(inOutR)
         isGhostL[i][0] =  isOutsideL and not isOutsideR
         isGhostR[i][0] = not isOutsideL and isOutsideR
      end

      -- loop over edges
      for i = localRange:lower(dir)-1, localRange:upper(dir)+1 do
         idxL[dir] = i
         idxR[dir] = i + 1
         if isGhostL[i][0] or isGhostR[i][0] then
            if (isGhostL[i][0]) then
               qOut:fill(indexer(idxL), qG)
               qOut:fill(indexer(idxR), qS)
            else -- right cell is ghost
               qOut:fill(indexer(idxR), qG)
               qOut:fill(indexer(idxL), qS)
            end
            for _, bc in ipairs(self._bcList) do
               bc(dir, tCurr, idxS, qS, qG) -- TODO: PASS COORDINATES
            end
         end
      end
   end

   self._isFirst = false
   return true, GKYL_MAX_DOUBLE
end

function StairSteppedBc:getDir()
   return self._dir
end

return StairSteppedBc
