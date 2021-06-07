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

         if ghostSkin or skinGhost then
            local idxGhost, idxSkin
            if (ghostSkin) then
               idxGhost, idxSkin = idxL, idxR
            else -- ghostSkin
               idxGhost, idxSkin = idxR, idxL
            end

            qOut:fill(indexer(idxGhost), qGhost)
            self._grid:setIndex(idxGhost)
            self._grid:cellCenter(xcGhost)

            qOut:fill(indexer(idxSkin), qSkin)
            self._grid:setIndex(idxSkin)
            self._grid:cellCenter(xcSkin)

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
