-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions at stair-stepped boundaries.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
local Alloc        = require "Lib.Alloc"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"

local isGhost = function (inOutPtr) return inOutPtr[1] < 0 end

local StairSteppedBc = Proto(UpdaterBase)

function StairSteppedBc:init(tbl)
   StairSteppedBc.super.init(self, tbl)
   local pfx = "Updater.StairSteppedBc: "

   self._grid = assert(tbl.onGrid, pfx.."Must specify 'onGrid'.")
   self._dir  = assert(tbl.dir,
                       pfx.."Must specify direction to apply BCs with 'dir'.")
   self._bcList = assert(tbl.boundaryConditions,
                         pfx.."Must specify 'boundaryConditions'.")
   self._inOut = assert(tbl.inOut, pfx.."Must specify mask field with 'inOut'.")

   self._perpRangeDecomp = LinearDecomp.LinearDecompRange {
      range      = self._grid:localRange():shorten(self._dir),
      numSplit   = self._grid:numSharedProcs(),
      threadComm = self:getSharedComm()
   }
end

function StairSteppedBc:_advance(tCurr, inFld, outFld)
   local qOut = assert(outFld[1],
                      "StairSteppedBc.advance: Must-specify an output field")
   local grid = self._grid
   local dir = self._dir
   local localRange = grid:localRange()

   -- Pointer and indices to be reused.
   local qGhost, qSkin = qOut:get(1), qOut:get(1)
   local idxL    = Lin.IntVec(grid:ndim())
   local idxR    = Lin.IntVec(grid:ndim())
   local indexer = qOut:genIndexer()

   local inOut          = self._inOut
   local inOutL, inOutR = inOut:get(1), inOut:get(1)
   local inOutIndexer      = inOut:genIndexer()

   local xcSkin   = Lin.Vec(self._grid:ndim())
   local xcGhost  = Lin.Vec(self._grid:ndim())

   -- Loop over directions orthogonal to 'dir'.
   local tId = self._grid:subGridSharedId()
   for idx in self._perpRangeDecomp:colMajorIter(tId) do
      idx:copyInto(idxL)
      idx:copyInto(idxR)

      -- Loop over 1D slice in `dir`.
      for i = localRange:lower(dir)-1, localRange:upper(dir)+1 do
         idxL[dir] = i
         inOut:fill(inOutIndexer(idxL), inOutL)
         local leftCellIsGhost = isGhost(inOutL)

         idxR[dir] = i + 1
         inOut:fill(inOutIndexer(idxR), inOutR)
         local rightCellIsGhost = isGhost(inOutR)

         local ghostSkin = leftCellIsGhost and (not rightCellIsGhost)
         local skinGhost = (not leftCellIsGhost) and rightCellIsGhost

         local idxSkin
         if ghostSkin or skinGhost then
            if (ghostSkin) then
               qOut:fill(indexer(idxL), qGhost)
               qOut:fill(indexer(idxR), qSkin)
               self._grid:setIndex(idxL);  self._grid:cellCenter(xcGhost);
               self._grid:setIndex(idxR);  self._grid:cellCenter(xcSkin);
               idxSkin = idxR
            else
               qOut:fill(indexer(idxR), qGhost)
               qOut:fill(indexer(idxL), qSkin)
               self._grid:setIndex(idxL);  self._grid:cellCenter(xcSkin);
               self._grid:setIndex(idxR);  self._grid:cellCenter(xcGhost);
               idxSkin = idxL
            end
            for _, bc in ipairs(self._bcList) do
               bc(dir, tCurr, idxSkin, qSkin, qGhost, xcGhost, xcSkin)
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
