-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the Spitzer collision frequency.
-- Eventually this updater may:
-- a) Return a cell-wise constant nu given a normalized initial normNu.
-- b) Return a nu expanded in the basis function given an initial normNu.
-- c) Build nu in SI units (either cell-wise constant or varying).
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase   = require "Updater.Base"
local LinearDecomp  = require "Lib.LinearDecomp"
local Proto         = require "Lib.Proto"
local SpitzerNuDecl = require "Updater.spitzerNuCalcData.SpitzerNuModDecl"
local xsys          = require "xsys"
local Lin           = require "Lib.Linalg"

-- Updater object.
local SpitzerCollisionality = Proto(UpdaterBase)

function SpitzerCollisionality:init(tbl)
   SpitzerCollisionality.super.init(self, tbl) -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SpitzerCollisionality: Must provide grid object using 'onGrid'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.SpitzerCollisionality: Must provide the configuration basis object using 'confBasis'.")

   self._cellConstNu = assert(
      tbl.useCellAverageNu, "Updater.SpitzerCollisionality: Must specify whether to use cell averaged collisionality in 'useCellAverageNu'.")

   self._isNormNu = tbl.willInputNormNu
   assert(self._isNormNu~=nil, "Updater.SpitzerCollisionality: Must specify whether normalized collisionality (normNu) will be inputed or not.")

   self._elemCharge = assert(
      tbl.elemCharge, "Updater.SpitzerCollisionality: Must specify elementary charge ('elemCharge') to build Spitzer collisionality.")
   self._epsilon0   = assert(
      tbl.epsilon0, "Updater.SpitzerCollisionality: Must specify vacuum permittivity ('epsilon0') to build Spitzer collisionality.")
   self._hBar       = assert(
      tbl.hBar, "Updater.SpitzerCollisionality: Must specify Planck's constant h divided by 2pi ('hBar') to build Spitzer collisionality.")
   local nuFracIn = tbl.nuFrac
   if nuFracIn then
      self._nuFrac = nuFracIn
   else
      self._nuFrac = 1.0
   end

   -- Dimension of configuration space.
   self._cDim = confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisC = confBasis:numBasis()

   if self._isNormNu then
      if self._cellConstNu then
        self._SpitzerNuCalc = SpitzerNuDecl.selectCellAvSpitzerNuScale(self._basisID, self._cDim, self._polyOrder)
      else
        self._SpitzerNuCalc = SpitzerNuDecl.selectSpitzerNuScale(self._basisID, self._cDim, self._polyOrder)
      end
   else
      if self._cellConstNu then
        self._SpitzerNuCalc = SpitzerNuDecl.selectCellAvSpitzerNuBuild(self._basisID, self._cDim, self._polyOrder)
      else
        self._SpitzerNuCalc = SpitzerNuDecl.selectSpitzerNuBuild(self._basisID, self._cDim, self._polyOrder)
      end
   end
   self.onGhosts = xsys.pickBool(false, tbl.onGhosts)

   self._BmagZero = Lin.Vec(self._numBasisC)
   for iC = 1,self._numBasisC do
      self._BmagZero[iC] = 0.0
   end
   
end

-- Advance method.
function SpitzerCollisionality:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local chargeSelf, massSelf    = inFld[1], inFld[2]
   local m0Self, vtSqSelf        = inFld[3], inFld[4]
   local chargeOther, massOther  = inFld[5], inFld[6]
   local m0Other, vtSqOther      = inFld[7], inFld[8]

   local normNu = inFld[9]
   local Bmag, BmagItr

   local firstInput, massFac = 0.0, 0.0

   local confIndexer   = m0Self:genIndexer()

   local m0SelfItr     = m0Self:get(1)
   local vtSqSelfItr   = vtSqSelf:get(1)
   local m0OtherItr    = m0Other:get(1)
   local vtSqOtherItr  = vtSqOther:get(1)

   local nuOut          = outFld[1]
   local nuOutItr       = nuOut:get(1)

   local confRange      = m0Self:localRange()
   if self.onGhosts then confRange = m0Self:localExtRange() end

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- Local thread ID.

   -- Fork logic here to avoid if-statement in space loop.
   if (inFld[10] ~= nil) then
      Bmag    = inFld[10]
      BmagItr = Bmag:get(1)
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)
 
         m0Self:fill(confIndexer(cIdx), m0SelfItr)
         vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         m0Other:fill(confIndexer(cIdx), m0OtherItr)
         vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)

         Bmag:fill(confIndexer(cIdx), BmagItr)
 
         nuOut:fill(confIndexer(cIdx), nuOutItr)
 
         self._SpitzerNuCalc(self._elemCharge, self._epsilon0, self._hBar, self._nuFrac, chargeSelf, massSelf, m0SelfItr:data(), vtSqSelfItr:data(), chargeOther, massOther, m0OtherItr:data(), vtSqOtherItr:data(), normNu, BmagItr:data(), nuOutItr:data())
      end
   else
      BmagItr = self._BmagZero
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)
 
         m0Self:fill(confIndexer(cIdx), m0SelfItr)
         vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         m0Other:fill(confIndexer(cIdx), m0OtherItr)
         vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)
 
         nuOut:fill(confIndexer(cIdx), nuOutItr)
 
         self._SpitzerNuCalc(self._elemCharge, self._epsilon0, self._hBar, self._nuFrac, chargeSelf, massSelf, m0SelfItr:data(), vtSqSelfItr:data(), chargeOther, massOther, m0OtherItr:data(), vtSqOtherItr:data(), normNu, BmagItr:data(), nuOutItr:data())
      end
   end
end

return SpitzerCollisionality
