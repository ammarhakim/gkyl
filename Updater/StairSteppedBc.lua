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

   local xcIn   = Lin.Vec(self._grid:ndim())
   local xcOut  = Lin.Vec(self._grid:ndim())

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
         local isGhostL =  isOutsideL and not isOutsideR
         local isGhostR = not isOutsideL and isOutsideR

         local idxIn
         if isGhostL or isGhostR then
            if (isGhostL) then
               qOut:fill(indexer(idxL), qG)
               qOut:fill(indexer(idxR), qS)
               self._grid:setIndex(idxL);  self._grid:cellCenter(xcOut);
               self._grid:setIndex(idxR);  self._grid:cellCenter(xcIn);
               idxIn = idxR
            else
               qOut:fill(indexer(idxR), qG)
               qOut:fill(indexer(idxL), qS)
               self._grid:setIndex(idxL);  self._grid:cellCenter(xcIn);
               self._grid:setIndex(idxR);  self._grid:cellCenter(xcOut);
               idxIn = idxL
            end
            for _, bc in ipairs(self._bcList) do
               bc(dir, tCurr, idxIn, qS, qG, xcOut, xcIn)
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
