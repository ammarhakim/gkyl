-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the primitive moments, u and vtSq=sqrt(T/m), given
-- the distribution function and its moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase     = require "Updater.Base"
local Lin             = require "Lib.Linalg"
local LinearDecomp    = require "Lib.LinearDecomp"
local Proto           = require "Lib.Proto"
local PrimMomentsDecl = require "Updater.primMomentsCalcData.PrimMomentsModDecl"
local xsys            = require "xsys"
local ffi = require "ffi"
local ffiC = ffi.C

ffi.cdef [[
  void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned c);
  void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out);
]]

-- function to check if operator option is correct
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" then
      return true
   end
   return false
end

-- Moments updater object.
local SelfPrimMoments = Proto(UpdaterBase)

function SelfPrimMoments:init(tbl)
   SelfPrimMoments.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SelfPrimMoments: Must provide grid object using 'onGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.SelfPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.SelfPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   local operator = assert(
      tbl.operator, "Updater.SelfPrimMoments: Must specify the collision operator (VmLBO or GkLBO for now) using 'operator'.")

   -- Dimension of phase space.
   self._pDim = phaseBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- Ensure sanity.
   assert(phaseBasis:polyOrder() == self._polyOrder,
          "Polynomial orders of phase and conf basis must match.")
   assert(phaseBasis:id() == self._basisID,
          "Type of phase and conf basis must match.")
   -- Determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions. Used to compute number of vector components.
   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   local uDim       = self._vDim  -- Dimensionality of flow velocity vector.
   self._kinSpecies = "Vm"        -- Vlasov-Maxwell species.
   self._isGkLBO    = false
   if isOperatorGood(operator) then
     if operator == "GkLBO" then
        self._isGkLBO    = true
        self._kinSpecies = "Gk"    -- Gyrokinetic species.
        uDim             = 1       -- A (vpar,mu) simulation has 3 physical velocity dimensions.
     end
   else
      assert(false, string.format(
                "SelfPrimMoments: Operator option must be 'VmLBO' or 'GkLBO'. Requested %s instead.", operator))
   end

   self._SelfPrimMomentsCalc = PrimMomentsDecl.selectSelfPrimMomentsCalc(self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
   
   self._binOpData = ffiC.new_binOpData_t(self._numBasisC*(uDim+1), 0) 
end

-- Advance method.
function SelfPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local uOut           = outFld[1]
   local vtSqOut        = outFld[2]

   local uOutItr        = uOut:get(1)
   local vtSqOutItr     = vtSqOut:get(1)

   -- Moments used for all polyOrders.
   local m0, m1 = inFld[1], inFld[2]
   local confIndexer   = m0:genIndexer()
   local m0Itr       = m0:get(1)
   local m1Itr       = m1:get(1)

   -- Boundary corrections.
   local cMomB, cEnergyB = inFld[4], inFld[5]
   local cMomBItr        = cMomB:get(1)
   local cEnergyBItr     = cEnergyB:get(1)

   local m2, m2Itr
   local m0Star, m0StarItr
   local m1Star, m1StarItr
   local m2Star, m2StarItr
   if self._polyOrder > 1 then
      m2         = inFld[3]
      m2Itr      = m2:get(1)
   else
      m0Star, m1Star, m2Star  = inFld[6], inFld[7], inFld[8]
      m0StarItr     = m0Star:get(1)
      m1StarItr     = m1Star:get(1)
      m2StarItr     = m2Star:get(1)
   end

   local confRange  = m0:localRange()
   if self.onGhosts then confRange = m0:localExtRange() end

   -- Construct ranges for nested loops.
   -- NOTE: Shared memory is only being used over configuration space.
   -- Similar to DistFuncMomentCalc, this is to avoid race conditions.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.

   -- polyOrder=1 and >1 each use separate velocity grid loops to
   -- avoid evaluating (if polyOrder==1) at each cell.
   if self._polyOrder > 1 then

      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)

         m0:fill(confIndexer(cIdx), m0Itr)
         m1:fill(confIndexer(cIdx), m1Itr)
         m2:fill(confIndexer(cIdx), m2Itr)
         cMomB:fill(confIndexer(cIdx), cMomBItr)
         cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)

         uOut:fill(confIndexer(cIdx), uOutItr)
         vtSqOut:fill(confIndexer(cIdx), vtSqOutItr)

         self._SelfPrimMomentsCalc(self._binOpData, m0Itr:data(), m1Itr:data(), m2Itr:data(), cMomBItr:data(), cEnergyBItr:data(), uOutItr:data(), vtSqOutItr:data())
      end

   else

      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)

         m0:fill(confIndexer(cIdx), m0Itr)
         m1:fill(confIndexer(cIdx), m1Itr)
         cMomB:fill(confIndexer(cIdx), cMomBItr)
         cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)
         m0Star:fill(confIndexer(cIdx), m0StarItr)
         m1Star:fill(confIndexer(cIdx), m1StarItr)
         m2Star:fill(confIndexer(cIdx), m2StarItr)

         uOut:fill(confIndexer(cIdx), uOutItr)
         vtSqOut:fill(confIndexer(cIdx), vtSqOutItr)

         self._SelfPrimMomentsCalc(self._binOpData, m0Itr:data(), m1Itr:data(), m0StarItr:data(), m1StarItr:data(), m2StarItr:data(), cMomBItr:data(), cEnergyBItr:data(), uOutItr:data(), vtSqOutItr:data())
      end

   end
end

return SelfPrimMoments
