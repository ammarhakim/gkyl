-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate BGK collisions for various species.
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local GaussQuadRules = require "Lib.GaussQuadRules"
local LinearDecomp   = require "Lib.LinearDecomp"
local Lin            = require "Lib.Linalg"
local Proto          = require "Proto"
local Range          = require "Lib.Range"
local Time           = require "Lib.Time"
local UpdaterBase    = require "Updater.Base"

-- System libraries.
local xsys = require "xsys"

-- BGK Collisions updater object.
local BgkCollisions = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function BgkCollisions:init(tbl)
   BgkCollisions.super.init(self, tbl) -- Setup base object.

   self._phaseGrid     = assert(tbl.onGrid,
      "Updater.BgkCollisions: Must provide phase space grid object using 'onGrid'")
   self._phaseBasis    = assert(tbl.phaseBasis,
      "Updater.BgkCollisions: Must provide phase space basis object using 'phaseBasis'")
   self._confGrid      = assert(tbl.confGrid,
      "Updater.BgkCollisions: Must provide configuration space grid object using 'confGrid'")
   self._confBasis     = assert(tbl.confBasis,
      "Updater.BgkCollisions: Must provide configuration space basis object using 'confBasis'")
   local varNuIn       = tbl.varyingNu           -- Specify if collisionality varies spatially.
   local cellConstNuIn = tbl.useCellAverageNu    -- Specify whether to use cell-wise constant collisionality.

  -- Dimension of phase space.
   self._pDim = self._phaseBasis:ndim()

   -- The default is spatially constant collisionality.
   if varNuIn==nil then
      self._varNu       = false
      self._cellConstNu = true
   else
      self._varNu       = varNuIn
      if cellConstNuIn == nil then
         self._cellConstNu = true
      else
         self._cellConstNu = cellConstNuIn
      end
   end

   if self._varNu then
      self._nuPtr, self._nuIdxr = nil, nil
   end

   local numConfDims   = self._confBasis:ndim()
   -- To obtain the cell average, multiply the zeroth coefficient by this factor.
   self._cellAvFac = 1.0/(math.sqrt(2.0^numConfDims))
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function BgkCollisions:_advance(tCurr, inFld, outFld)
   local grid = self._phaseGrid
   local numPhaseBasis = self._phaseBasis:numBasis()
   local pDim = self._pDim
   -- Get the inputs and outputs.
   local fIn           = assert(inFld[1],
      "BgkCollisions.advance: Must specify an input distribution function field as input[1]")
   local sumNufMaxwell = assert(inFld[2],
      "BgkCollisions.advance: Must specify sum(nu*Maxwellian) field as input[2]")
   local sumNuIn       = assert(inFld[3],
      "BgkCollisions.advance: Must specify the sum of collisionalities as input[3]")
   local nuFrac = 1.0
   if inFld[4] then
      nuFrac = inFld[4]
   end

   local fRhsOut = assert(outFld[1], "BgkCollisions.advance: Must specify an output field")

   local fInItr           = fIn:get(1)
   local sumNufMaxwellItr = sumNufMaxwell:get(1)

   local sumNu = 0.0    -- Assigned below if cellConstNu=true.
   if self._varNu then
      self._sumNuPtr  = sumNuIn:get(1)
      self._sumNuIdxr = sumNuIn:genIndexer()
   else
      sumNu = sumNuIn
   end

   local fRhsOutItr   = fRhsOut:get(1)
   local phaseRange = fRhsOut:localRange()
   -- Get the range to loop over the domain.
   local phaseIndexer = fRhsOut:genIndexer()

   -- Get the interface for setting global CFL frequencies
   local cflRateByCell = self._cflRateByCell
   local cflRateByCellIdxr = cflRateByCell:genIndexer()
   local cflRateByCellPtr = cflRateByCell:get(1)

   -- Construct range for shared memory.
   local phaseRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(pDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.

   -- Phase space loop.
   for pIdx in phaseRangeDecomp:rowMajorIter(tId) do
      fIn:fill(phaseIndexer(pIdx), fInItr)
      sumNufMaxwell:fill(phaseIndexer(pIdx), sumNufMaxwellItr)
      fRhsOut:fill(phaseIndexer(pIdx), fRhsOutItr)

      cflRateByCell:fill(cflRateByCellIdxr(pIdx), cflRateByCellPtr)

      if self._cellConstNu then
         -- This code assumes nu is cell-wise constant.
         if self._varNu then
            sumNuIn:fill(self._sumNuIdxr(pIdx), self._sumNuPtr)    -- Get pointer to sumNu field.
            sumNu = self._sumNuPtr[1]*self._cellAvFac
         end
         for k = 1, numPhaseBasis do
            if sumNufMaxwellItr[k] == sumNufMaxwellItr[k] then -- NaN check.
               fRhsOutItr[k] = fRhsOutItr[k] + nuFrac*(sumNufMaxwellItr[k] - sumNu*fInItr[k])
            end
         end
	 cflRateByCellPtr:data()[0] = cflRateByCellPtr:data()[0] + sumNu
      end
   end
end

return BgkCollisions
