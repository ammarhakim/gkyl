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

local isGhost = function (inOutPtr)
   return inOutPtr[1] < 0
end

function StairSteppedBc:init(tbl)
   StairSteppedBc.super.init(self, tbl) -- setup base object

   self._grid = assert(tbl.onGrid, "Updater.StairSteppedBc: Must specify grid to use with 'onGrid''")
   self._dir  = assert(tbl.dir, "Updater.StairSteppedBc: Must specify direction to apply BCs with 'dir'")

   self._bcList = assert(
   tbl.boundaryConditions, "Updater.StairSteppedBc: Must specify boundary conditions to apply with 'boundaryConditions'")

   self._localRange = tbl.onGrid:localRange()
   local localRange = self._localRange

   self._inOut = assert(tbl.inOut, "Updater.StairSteppedBc: Must specify mask field with 'inOut'")

   self._perpRangeDecomp = LinearDecomp.LinearDecompRange {
      range      = localRange:shorten(self._dir), -- range orthogonal to 'dir'
      numSplit   = self._grid:numSharedProcs(),
      threadComm = self:getSharedComm()
   }
end

function StairSteppedBc:_advance(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "StairSteppedBc.advance: Must-specify an output field")
   
   local dir = self._dir
   local localRange = grid:localRange()

   -- Pointer and indices to be (re)used.
   local qGhost, qSkin = qOut:get(1), qOut:get(1)
   local idxL    = Lin.IntVec(grid:ndim())
   local idxR    = Lin.IntVec(grid:ndim())
   local indexer = qOut:genIndexer()

   local inOut          = self._inOut
   local inOutL, inOutR = inOut:get(1), inOut:get(1)
   local inOutIdxr      = inOut:genIndexer()

   local xcIn   = Lin.Vec(self._grid:ndim())
   local xcOut  = Lin.Vec(self._grid:ndim())

   -- Loop over directions orthogonal to 'dir'.
   local tId = self._grid:subGridSharedId()
   for idx in self._perpRangeDecomp:colMajorIter(tId) do
      idx:copyInto(idxL)
      idx:copyInto(idxR)

      -- Loop over 1D slice in `dir`.
      for i = localRange:lower(dir)-1, localRange:upper(dir)+1 do
         idxL[dir] = i
         inOut:fill(inOutIdxr(idxL), inOutL)
         local leftCellIsGhost = isGhost(inOutL)

         idxR[dir] = i + 1
         inOut:fill(inOutIdxr(idxR), inOutR)
         local rightCellIsGhost = isGhost(inOutR)

         local ghostSkin = leftCellIsGhost and (not rightCellIsGhost)
         local skinGhost = (not leftCellIsGhost) and rightCellIsGhost

         local idxIn
         if ghostSkin or skinGhost then
            if (ghostSkin) then
               qOut:fill(indexer(idxL), qGhost)
               qOut:fill(indexer(idxR), qSkin)
               self._grid:setIndex(idxL);  self._grid:cellCenter(xcOut);
               self._grid:setIndex(idxR);  self._grid:cellCenter(xcIn);
               idxIn = idxR
            else
               qOut:fill(indexer(idxR), qGhost)
               qOut:fill(indexer(idxL), qSkin)
               self._grid:setIndex(idxL);  self._grid:cellCenter(xcIn);
               self._grid:setIndex(idxR);  self._grid:cellCenter(xcOut);
               idxIn = idxL
            end
            for _, bc in ipairs(self._bcList) do
               bc(dir, tCurr, idxIn, qSkin, qGhost, xcOut, xcIn)
            end
         end
      end
   end

   return true, GKYL_MAX_DOUBLE
end

function StairSteppedBc:getDir()
   return self._dir
end

return StairSteppedBc
