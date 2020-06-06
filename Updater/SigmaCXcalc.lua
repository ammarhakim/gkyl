-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate SigmaCX for the Pauls CX model
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase        = require "Updater.Base"
local LinearDecomp       = require "Lib.LinearDecomp"
local Proto              = require "Lib.Proto"
local ChargeExchangeDecl = require "Updater.chargeExchangeCalcData.ChargeExchangeModDecl"
local xsys               = require "xsys"
local Lin                = require "Lib.Linalg"
local Time               = require "Lib.Time"

-- Charge exchange collisions updater object.
local SigmaCX = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function SigmaCX:init(tbl)
   SigmaCX.super.init(self, tbl) -- setup base object

   self._onGrid     = assert(tbl.onGrid,
			     "Updater.SigmaCX: Must provide grid object using 'onGrid'")
   self._confBasis  = assert(tbl.confBasis,
			     "Updater.SigmaCX: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.SigmaCX: Must provide phase space basis object using 'phaseBasis'")
   self._kineticSpecies = assert(tbl.kineticSpecies,
			   "Updater.SigmaCX: Must provide solver type (Vm or Gk) using 'kineticSpecies'")
   self._a = assert(tbl.a,
		    "Updater.SigmaCX: Must provide fitting constant a using 'a'")
   self._b = assert(tbl.b,
		    "Updater.SigmaCX: Must provide fitting constant b using 'b'")
   
   -- Dimension of spaces.
   self._cDim = self._confBasis:ndim()
   self._pDim = self._phaseBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Basis name and polynomial order.
   self._basisID   = self._phaseBasis:id()
   self._polyOrder = self._phaseBasis:polyOrder()

   -- Number of basis functions.
   self._numBasis = self._confBasis:numBasis()

   -- Define CX cross section calculation
   if self._kineticSpecies == "Vm" then
      self._calcSigmaCX = ChargeExchangeDecl.VmSigmaCX(self._basisID, self._cDim, self._vDim, self._polyOrder)
   elseif self._kineticSpecies == "Gk" then 
      self._calcSigmaCX = ChargeExchangeDecl.GkSigmaCX(self._basisID, self._cDim, self._vDim, self._polyOrder)
   else
      print("Updater.SigmaCX: 'kineticSpecies must be 'Vm' or 'Gk'")
   end
   
   self.onGhosts = xsys.pickBool(false, tbl.onGhosts)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function SigmaCX:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local numPhaseBasis = self._phaseBasis:numBasis()
   local cDim = self._cDim

   local uIon       = assert(inFld[1], "SigmaCX.advance: Must specify ion fluid velocity as input[1]")
   local uNeut      = assert(inFld[2], "SigmaCX.advance: Must specify neutral fluid velocity as input[2]")
   local vtSqIon    = assert(inFld[3], "SigmaCX.advance: Must specify ion squared thermal velocity as input[3]")
   local vtSqNeut   = assert(inFld[4], "SigmaCX.advance: Must specify neutral squared thermal velocity as input[4]")
   local sigmaCX    = assert(outFld[1], "SigmaCX.advance: Must specify an output field")
   
   local confIndexer  = uIon:genIndexer()

   local uIonItr      = uIon:get(1)
   local uNeutItr     = uNeut:get(1)
   local vtSqIonItr   = vtSqIon:get(1)
   local vtSqNeutItr  = vtSqNeut:get(1)
 
   local sigmaCXItr   = sigmaCX:get(1)

   local confRange = vtSqIon:localRange()

   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.
   
   -- Phase space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(cIdx)
      
      uIon:fill(confIndexer(cIdx), uIonItr)
      uNeut:fill(confIndexer(cIdx), uNeutItr)      
      vtSqIon:fill(confIndexer(cIdx), vtSqIonItr)
      vtSqNeut:fill(confIndexer(cIdx), vtSqNeutItr)
      sigmaCX:fill(confIndexer(cIdx), sigmaCXItr)

      self._calcSigmaCX(self._a, self._b, uIonItr:data(), uNeutItr:data(), vtSqIonItr:data(), vtSqNeutItr:data(), sigmaCXItr:data())
     
   end
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function SigmaCX:evalMomTime() return self._tmEvalMom end

return SigmaCX
