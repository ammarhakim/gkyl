-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate relative velocity v^*,
-- and return product of v^*, density and distF
-- for the Pauls CX model
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local Proto = require "Lib.Proto"
local RelativeVelocityDecl = require "Updater.chargeExchangeCalcData.RelativeVelocityModDecl"
local xsys = require "xsys"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

-- Charge exchange collisions updater object.
local VrelProductCX = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function VrelProductCX:init(tbl)
   VrelProductCX.super.init(self, tbl) -- setup base object

   self._onGrid     = assert(tbl.onGrid,
			     "Updater.VrelProductCX: Must provide grid object using 'onGrid'")
   self._confBasis  = assert(tbl.confBasis,
			     "Updater.VrelProductCX: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.VrelProductCX: Must provide phase space basis object using 'phaseBasis'")
   self._kineticSpecies = assert(tbl.kineticSpecies,
			   "Updater.VrelProductCX: Must provide solver type (Vm or Gk) using 'kineticSpecies'")   
   -- Dimension of spaces.
   self._pDim = self._phaseBasis:ndim()
   self._cDim = self._confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Basis name and polynomial order.
   self._basisID   = self._phaseBasis:id()
   self._polyOrder = self._phaseBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisP = self._phaseBasis:numBasis()

   -- Define cell center
   self.xc  = Lin.Vec(self._pDim)

   -- Define relative velocity calculation
   if self._kineticSpecies == "Vm" then
      self._calcVrelProductCX = RelativeVelocityDecl.VmVrelProdCX(self._basisID, self._cDim, self._vDim, self._polyOrder)
   elseif self._kineticSpecies == "Gk" then 
      self._calcVrelProductCX = RelativeVelocityDecl.GkVrelProdCX(self._basisID, self._cDim, self._vDim, self._polyOrder)
   else
      print("Updater.SigmaCX: 'kineticSpecies must be 'Vm' or 'Gk'")
   end
   
   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function VrelProductCX:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local pDim = self._pDim

   local m0     = assert(inFld[1], "VrelProductCX.advance: Must specify particle density as input[1]")
   local u      = assert(inFld[2], "VrelProductCX.advance: Must specify fluid velocity as input[2]")
   local vtSq   = assert(inFld[3], "VrelProductCX.advance: Must specify squared thermal velocity as input[3]")
   local fOther = assert(inFld[4], "VrelProductCX.advance: Must specify distF of other species as input[4]")
   local prodCX = assert(outFld[1], "VrelProductCX.advance: Must specify an output field")
   
   local confIndexer  = vtSq:genIndexer()
   local phaseIndexer = prodCX:genIndexer()

   local m0Itr     = m0:get(1)   
   local uItr      = u:get(1)
   local vtSqItr   = vtSq:get(1)
   local fOtherItr = fOther:get(1)   
   local prodCXItr = prodCX:get(1)

   local phaseRange = prodCX:localRange()

   -- Phase space loop 
   for pIdx in phaseRange:rowMajorIter() do
      grid:setIndex(pIdx)
      grid:cellCenter(self.xc)

      m0:fill(confIndexer(pIdx), m0Itr)
      u:fill(confIndexer(pIdx), uItr)
      vtSq:fill(confIndexer(pIdx), vtSqItr)
      fOther:fill(phaseIndexer(pIdx), fOtherItr)
      prodCX:fill(phaseIndexer(pIdx),prodCXItr)

      self._calcVrelProductCX(self.xc:data(), m0Itr:data(), uItr:data(), vtSqItr:data(), fOtherItr:data(), prodCXItr:data())
     
   end
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function VrelProductCX:evalMomTime() return self._tmEvalMom end

return VrelProductCX
