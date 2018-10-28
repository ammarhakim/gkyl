-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the Spitzer collision frequency.
-- Eventually this updater may:
-- a) Return a cell-wise constant nu given a normalized initial normNu.
-- b) Return a nu expanded in the basis function given an initial normNi.
-- c) Construct nu in SI units (either cell-wise constant or varying).
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase   = require "Updater.Base"
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

   self._normNu      = assert(tbl.normalizedNu, "Updater.SpitzerCollisionality: Use 'normalizedNu' to specify the collisionality normalized by (T_0^(3/2)/n_0), where these are values somewhere in the domain at t=0.")
   self._cellConstNu = assert(
      tbl.useCellAverageNu, "Updater.SpitzerCollisionality: Must specify whether to use cell averaged collisionality in 'useCellAverageNu'.")
   local massIn      = assert(tbl.mass, "Updater.SpitzerCollisionality: Must provide the species mass with 'mass'.")

   self._rMassR3d2   = 1.0/math.sqrt(massIn^3)

   -- Dimension of configuration space.
   self._cDim = confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisC = confBasis:numBasis()

   if self._cellConstNu then
     self._SpitzerNuCalc = SpitzerNuDecl.selectCellAvSpitzerNuCalc(self._basisID, self._cDim, self._polyOrder)
   else
     self._SpitzerNuCalc = SpitzerNuDecl.selectSpitzerNuCalc(self._basisID, self._cDim, self._polyOrder)
   end
   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
   
end

-- Advance method.
function SpitzerCollisionality:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid

   local nuOut          = outFld[1]

   local nuOutIndexer   = nuOut:genIndexer()
   local nuOutItr       = nuOut:get(1)

   local m0Fld, vtSqFld = inFld[1], inFld[2]

   local m0FldIndexer   = m0Fld:genIndexer()
   local m0FldItr       = m0Fld:get(1)
   local vtSqFldIndexer = vtSqFld:genIndexer()
   local vtSqFldItr     = vtSqFld:get(1)

   local confRange      = m0Fld:localRange()
   if self.onGhosts then confRange = m0Fld:localExtRange() end

   for confIdx in confRange:colMajorIter() do
      grid:setIndex(confIdx)

      m0Fld:fill(m0FldIndexer(confIdx), m0FldItr)
      vtSqFld:fill(vtSqFldIndexer(confIdx), vtSqFldItr)
      nuOut:fill(nuOutIndexer(confIdx), nuOutItr)

      self._SpitzerNuCalc(self._normNu, self._rMassR3d2, m0FldItr:data(), vtSqFldItr:data(), nuOutItr:data())
   end

   return true, GKYL_MAX_DOUBLE
end

return SpitzerCollisionality
