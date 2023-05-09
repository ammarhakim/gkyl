-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate SigmaCX for the Pauls CX model
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local ChargeExchangeDecl = require "Updater.chargeExchangeCalcData.ChargeExchangeModDecl"
local xsys = require "xsys"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

-- Charge exchange collisions updater object.
local ChargeExchange = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function ChargeExchange:init(tbl)
   ChargeExchange.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid,
			     "Updater.SigmaCX: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis,
			     "Updater.SigmaCX: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.SigmaCX: Must provide phase space basis object using 'phaseBasis'")
   -- self._kineticSpecies = assert(tbl.kineticSpecies,
   -- 			   "Updater.SigmaCX: Must provide solver type (Vm or Gk) using 'kineticSpecies'")
   self._a = assert(tbl.a,
		    "Updater.SigmaCX: Must provide fitting constant a using 'a'")
   self._b = assert(tbl.b,
		    "Updater.SigmaCX: Must provide fitting constant b using 'b'")
   self._charge = tbl.charge
   
   -- Dimension of spaces.
   self._cDim = self._confBasis:ndim()
   self._pDim = self._phaseBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Basis name and polynomial order.
   self._basisID = self._phaseBasis:id()
   self._polyOrder = self._phaseBasis:polyOrder()
   self.idxP = Lin.IntVec(self._pDim)
   
   -- Number of basis functions.
   self._numBasis = self._confBasis:numBasis()

   -- Define CX cross section calculation
   if self._charge == 0 then
      self._calcSigmaCX = ChargeExchangeDecl.SigmaCX(self._basisID, self._cDim, self._vDim, self._polyOrder)
   end
   
   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function ChargeExchange:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local numConfBasis = self._confBasis:numBasis()
   local numPhaseBasis = self._phaseBasis:numBasis()
   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local m0          = assert(inFld[1], "SigmaCX.advance: Must specify neutral particle density as input[1]")
   local uIon        = assert(inFld[2], "SigmaCX.advance: Must specify ion fluid velocity as input[2]")
   local uNeut       = assert(inFld[3], "SigmaCX.advance: Must specify neutral fluid velocity as input[3]")
   local vtSqIon     = assert(inFld[4], "SigmaCX.advance: Must specify ion squared thermal velocity as input[4]")
   local vtSqIonMin  = assert(inFld[5], "SigmaCX.advance: Must specify ion squared thermal velocity as input[5]")
   local vtSqNeut    = assert(inFld[6], "SigmaCX.advance: Must specify neutral squared thermal velocity as input[6]")
   local vtSqNeutMin = assert(inFld[7], "SigmaCX.advance: Must specify neutral squared thermal velocity as input[7]")
   local vSigmaCX    = assert(outFld[1], "SigmaCX.advance: Must specify an output field")
   
   local confIndexer  = vtSqIon:genIndexer()

   local uIonDim = uIon:numComponents()/numConfBasis -- number of dimensions in uIon, needed for GkSpecies

   local m0Itr = m0:get(1)
   local uIonItr = uNeut:get(1) -- use shape of uNeut
   local uParIonItr = uIon:get(1)
   local uNeutItr = uNeut:get(1)
   local vtSqIonItr = vtSqIon:get(1)
   local vtSqNeutItr = vtSqNeut:get(1)
   local vSigmaCXItr = vSigmaCX:get(1)
   
   local confRange = vtSqIon:localRange()

   -- Get the interface for setting global CFL frequencies
   local cflRateByCell     = self._cflRateByCell
   local cflRateByCellIdxr = cflRateByCell:genIndexer()
   local cflRateByCellPtr  = cflRateByCell:get(1)
   local cflRange = cflRateByCell:localRange()
   local velRange = cflRange:selectLast(vDim)

   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.

   local cflRate
   -- Conf space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(cIdx)

      if uIonDim < vDim and uIonDim == 1 and cDim == 1 then
	 uIon:fill(confIndexer(cIdx), uParIonItr)
	 for k = 1, numConfBasis do
	    uIonItr[k] = uParIonItr[k]
	 end
	 for d = 2, vDim do
	    for k = 1, numConfBasis do
	       uIonItr[numConfBasis*(d-1)+k] = 0.0
	    end
	 end
      elseif uIonDim < vDim and uIonDim == 1 and cDim == 3 then
	 uIon:fill(confIndexer(cIdx), uParIonItr)
	 for d = 1, vDim-1 do
	    for k = 1, numConfBasis do
	       uIonItr[numConfBasis*(d-1)+k] = 0.0
	    end
	 end
	 for k = 1, numConfBasis do
	    uIonItr[numConfBasis*(vDim-1)+k] = uParIonItr[k]
	 end
      elseif uIonDim == vDim then
	 uIon:fill(confIndexer(cIdx), uIonItr)
      else
	 print("Updater.SigmaCX: incorrect uIonDim")
      end
      m0:fill(confIndexer(cIdx), m0Itr)
      uNeut:fill(confIndexer(cIdx), uNeutItr)
      vtSqIon:fill(confIndexer(cIdx), vtSqIonItr)
      vtSqNeut:fill(confIndexer(cIdx), vtSqNeutItr)
      vSigmaCX:fill(confIndexer(cIdx), vSigmaCXItr)

      cflRate = self._calcSigmaCX(self._a, self._b, m0Itr:data(), uIonItr:data(), uNeutItr:data(), vtSqIonItr:data(), vtSqIonMin, vtSqNeutItr:data(), vtSqNeutMin, vSigmaCXItr:data())

      for vIdx in velRange:rowMajorIter() do
      	 -- Construct the phase space index ot of the configuration
      	 -- space and velocity space indices
         cIdx:copyInto(self.idxP)
         for d = 1, vDim do self.idxP[self._cDim+d] = vIdx[d] end
      	 cflRateByCell:fill(cflRateByCellIdxr(self.idxP), cflRateByCellPtr)
      	 cflRateByCellPtr:data()[0] = cflRateByCellPtr:data()[0] + cflRate
      end
     
   end
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function ChargeExchange:evalMomTime() return self._tmEvalMom end

return ChargeExchange
