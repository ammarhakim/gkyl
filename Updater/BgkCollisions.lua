-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate BGK collisions for various species.
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local GaussQuadRules = require "Lib.GaussQuadRules"
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

   -- The default is spatially constant collisionality.
   if varNuIn==nil then
      self._varNu       = false
      self._cellConstNu = true
   else
      self._varNu       = varNuIn
      if cellConstNuIn==nil then
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
   local numPhaseBasis = self._phaseBasis:numBasis()

   -- Get the inputs and outputs.
   local fIn      = assert(inFld[1],
      "BgkCollisions.advance: Must specify an input distribution function field as input[1]")
   local fMaxwell = assert(inFld[2],
      "BgkCollisions.advance: Must specify the Maxwellian distribution function field as input[2]")
   local nuIn     = assert(inFld[3],
      "BgkCollisions.advance: Must specify the collisionality as input[3]")
   local nuFrac = 1.0
   if inFld[4] then
      nuFrac = inFld[4]
   end

   local fRhsOut = assert(outFld[1], "BgkCollisions.advance: Must specify an output field")

   local fInItr      = fIn:get(1)
   local fMaxwellItr = fMaxwell:get(1)

   local nu = 0.0    -- Assigned below if cellConstNu=true.
   if self._varNu then
      self._nuPtr  = nuIn:get(1)
      self._nuIdxr = nuIn:genIndexer()
   else
      nu = nuIn
   end

   local fRhsOutItr  = fRhsOut:get(1)

   -- Get the range to loop over the domain.
   local phaseIndexer = fRhsOut:genIndexer()

   -- Phase space loop.
   for phaseIdx in fRhsOut:localRangeIter() do
      fIn:fill(phaseIndexer(phaseIdx), fInItr)
      fMaxwell:fill(phaseIndexer(phaseIdx), fMaxwellItr)
      fRhsOut:fill(phaseIndexer(phaseIdx), fRhsOutItr)

      if self._cellConstNu then
         -- This code assumes nu is cell-wise constant.
         if self._varNu then
            nuIn:fill(self._nuIdxr(phaseIdx), self._nuPtr)    -- Get pointer to nu field.
            nu = self._nuPtr[1]*self._cellAvFac
         end
         for k = 1, numPhaseBasis do
            if fMaxwellItr[k] == fMaxwellItr[k] then -- NaN check.
               fRhsOutItr[k] = fRhsOutItr[k] + nuFrac*nu*(fMaxwellItr[k] - fInItr[k])
            end
         end
      end
   end
end

return BgkCollisions
