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
   local epsilon0   = assert(
      tbl.epsilon0, "Updater.SpitzerCollisionality: Must specify vacuum permittivity ('epsilon0') to build Spitzer collisionality.")
   self._epsilon0Sq = epsilon0^2

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
   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
   
end

-- Advance method.
function SpitzerCollisionality:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local mass, charge   = inFld[1], inFld[2]
   local m0Fld, vtSqFld = inFld[3], inFld[4]

   local massSq   = mass^2
   local chargeR4 = charge^4

   local firstInput, massFac = 0.0, 0.0
   if self._isNormNu then
      firstInput = inFld[5]    -- Collisionality normalized by T_0^(3/2)/n_0.
      massFac    = 1.0/math.sqrt(mass^3)
   else
      firstInput = self._elemCharge    -- Elementary charge value.
      massFac    = mass
   end

   local confIndexer   = m0Fld:genIndexer()
   local m0FldItr       = m0Fld:get(1)
   local vtSqFldItr     = vtSqFld:get(1)

   local nuOut          = outFld[1]
   local nuOutItr       = nuOut:get(1)

   local confRange      = m0Fld:localRange()
   if self.onGhosts then confRange = m0Fld:localExtRange() end

   -- construct ranges for nested loops
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- local thread ID

   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(cIdx)

      m0Fld:fill(confIndexer(cIdx), m0FldItr)
      vtSqFld:fill(confIndexer(cIdx), vtSqFldItr)
      nuOut:fill(confIndexer(cIdx), nuOutItr)

      self._SpitzerNuCalc(firstInput, massFac, massSq, chargeR4, self._epsilon0Sq, m0FldItr:data(), vtSqFldItr:data(), nuOutItr:data())
   end
end

return SpitzerCollisionality
