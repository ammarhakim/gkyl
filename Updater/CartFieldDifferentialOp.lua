-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the result of a differential operator
-- acting on a Cartesian field quantity.
-- At the moment we are interested in the second derivative in 1D.
--
-- One could envision this updater computing 1st, 2nd derivatives, or Laplacians
-- in arbitrary directions and dimensions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase        = require "Updater.Base"
local Proto              = require "Lib.Proto"
local DifferentialOpDecl = require "Updater.differentialOpCalcData.CartFieldDifferentialOpModDecl"
local xsys               = require "xsys"
local Lin                = require "Lib.Linalg"
local LinearDecomp       = require "Lib.LinearDecomp"

-- Function to check if operation requested is correct/supported.
local function isOpNameGood(nm)
   if nm == "Dxx" then
      return true
   end
   return false
end

-- Updater object.
local CartFieldDifferentialOp = Proto(UpdaterBase)

function CartFieldDifferentialOp:init(tbl)
   CartFieldDifferentialOp.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.CartFieldDifferentialOp: Must provide grid object using 'onGrid'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.CartFieldDifferentialOp: Must provide the configuration basis object using 'confBasis'.")

   local op = assert(
      tbl.operation, "Updater.CartFieldDifferentialOp: Must provide an operation using 'operation'.")

   -- Dimension of configuration space.
   self._cDim = confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisC = confBasis:numBasis()

   -- By default, update all directions.
   self._updateDirs = {}
   for d = 1, self._cDim do
      self._updateDirs[d] = d
   end
   -- Read in which directions we are to update.
   if tbl.updateDirections then
      self._updateDirs = tbl.updateDirections
   end

   -- Don't perform surface updates in these directions (analogous to zeroFluxFlags in HyperDisCont).
   self._zeroSurfUpdateFlags = {}
   for d = 1, self._cDim do
      self._zeroSurfUpdateFlags[d] = false
   end
   local anyZeroSurfUpdateFlagsTrue = false
   if tbl.zeroSurfUpdateDirections then
      for i, d in ipairs(tbl.zeroSurfUpdateDirections) do
         self._zeroSurfUpdateFlags[d] = true
         anyZeroSurfUpdateFlagsTrue = true
      end
   end

   if isOpNameGood(op) then
      self._OperatorCalcSurf         = DifferentialOpDecl.selectOperatorSurf(op,self._basisID, self._cDim, self._polyOrder)
      self._OperatorCalcVol          = DifferentialOpDecl.selectOperatorVol(op,self._basisID, self._cDim, self._polyOrder)
      if anyZeroSurfUpdateFlagsTrue then
         self._OperatorCalcBoundarySurf = DifferentialOpDecl.selectOperatoruBoundarySurf(op,self._basisID, self._cDim, self._polyOrder)
      end
   else
      assert(false, string.format(
                "CartFieldDifferentialOp: Invalid differential operator. Requested %s instead.", op))
   end
   
   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   self._isFirst         = true
   self._perpRangeDecomp = {}    -- Perp ranges in each direction.
   
end

-- Advance method.
function CartFieldDifferentialOp:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local fIn        = inFld[1]
   local fInIdxr    = fIn:genIndexer()
   local fInM, fInP = fIn:get(1), fIn:get(1)

   local dfOut          = outFld[1]
   local dfOutIdxr      = dfOut:genIndexer()
   local dfOutM, dfOutP = dfOut:get(1), dfOut:get(1)

   -- Currently assume 1D:
   local ndim = 1
   local idxp, idxm = Lin.IntVec(ndim), Lin.IntVec(ndim)    -- Index on right/left.
   local dxp, dxm   = Lin.Vec(ndim), Lin.Vec(ndim)          -- Cell shape on right/left.
   local xcp, xcm   = Lin.Vec(ndim), Lin.Vec(ndim)          -- Cell center on right/left.

   local localRange                 = fIn:localRange()
   if self.onGhosts then localRange = fIn:localExtRange() end

   local tId = grid:subGridSharedId()    -- Local thread ID.

   local firstDir = true

   -- Iterate through updateDirs backwards so that a zeroSurfUpdate dir is first in kinetics.
   for i = #self._updateDirs, 1, -1 do
      local dir = self._updateDirs[i]

      -- Lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell).
      local dirLoIdx, dirUpIdx         = localRange:lower(dir), localRange:upper(dir)+1
      local dirLoSurfIdx, dirUpSurfIdx = dirLoIdx, dirUpIdx

      -- Compute loop bounds for directions without a surface update.
      if self._zeroSurfUpdateFlags[dir] then
         local dirGlobalLoIdx, dirGlobalUpIdx = globalRange:lower(dir), globalRange:upper(dir)+1
         if dirLoIdx == dirGlobalLoIdx then
            dirLoSurfIdx = dirLoIdx+1
         end
         if dirUpIdx == dirGlobalUpIdx then
            dirUpSurfIdx = dirUpIdx-1
         end
      end

      if self._isFirst then
         self._perpRangeDecomp[dir] = LinearDecomp.LinearDecompRange {
            range      = localRange:shorten(dir),    -- Range orthogonal to 'dir'.
            numSplit   = grid:numSharedProcs(),
            threadComm = self:getSharedComm()
         }
      end
      local perpRangeDecomp = self._perpRangeDecomp[dir]

      for idx in perpRangeDecomp:colMajorIter(tId) do
         idx:copyInto(idxp); idx:copyInto(idxm)

         for i = dirLoIdx, dirUpIdx do    -- This loop is over edges.
            idxm[dir], idxp[dir] = i-1, i -- Cell left/right of edge 'i'.
   
            grid:setIndex(idxm)
            for d = 1, ndim do dxm[d] = grid:dx(d) end
            grid:cellCenter(xcm)
       
            grid:setIndex(idxp)
            for d = 1, ndim do dxp[d] = grid:dx(d) end
            grid:cellCenter(xcp)
       
            fIn:fill(fInIdxr(idxm), fInM)
            fIn:fill(fInIdxr(idxp), fInP)
      
            dfOut:fill(dfOutIdxr(idxm), dfOutM)
            dfOut:fill(dfOutIdxr(idxp), dfOutP)
   
            if firstDir and i<=dirUpIdx-1 then
               self._OperatorCalcVol(xcp:data(), dxp:data(), idxp:data(), fInP:data(), dfOutP:data())
            end
            if i >= dirLoSurfIdx and i <= dirUpSurfIdx then
               self._OperatorCalcSurf[dir](
                  xcm:data(), xcp:data(), dxm:data(), dxp:data(), idxm:data(), idxp:data(), fInM:data(), fInP:data(), dfOutM:data(), dfOutP:data())
            else
               if self._zeroSurfUpdateFlags[dir] then
                  -- We need to give operators a chance to apply partial
                  -- surface updates even when the zeroFlux BCs have been
                  -- applied.
                  self._OperatorCalcBoundarySurf[dir](
                     xcm:data(), xcp:data(), dxm:data(), dxp:data(), idxm:data(), idxp:data(), fInM:data(), dInP:data(), dfOutM:data(), dfOutP:data())
               end
            end
         end
      end
      firstDir = false
   end

   self._isFirst = false
end

return CartFieldDifferentialOp
