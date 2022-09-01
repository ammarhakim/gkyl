-- Gkyl ------------------------------------------------------------------------
--
-- Updater for Voronov ionization calculations:
--     temperature (vtSqIz) for electrons or Voronov reaction rate (coefIz)
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase    = require "Updater.Base"
local Proto          = require "Lib.Proto"
local IonizationDecl = require "Updater.ionizationCalcData.IonizationModDecl"
local xsys = require "xsys"
local Lin  = require "Lib.Linalg"
local Time = require "Lib.Time"

-- Voronov Collisions updater object.
local Ionization = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function Ionization:init(tbl)
   Ionization.super.init(self, tbl) -- setup base object
   
   self._onGrid = assert(tbl.onGrid,
			  "Updater.Ionization: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis,
			    "Updater.Ionization: Must provide configuration space basis object using 'confBasis'")
   self._phaseGrid = assert(tbl.phaseGrid,
			  "Updater.Ionization: Must provide grid object using 'phaseGrid'")
   self._phaseBasis = assert(tbl.phaseBasis,
			    "Updater.Ionization: Must provide configuration space basis object using 'phaseBasis'")
   self._elcMass = assert(tbl.elcMass,
			  "Updater.Ionization: Must provide electron mass using 'elcMass'")
   self._elemCharge = assert(tbl.elemCharge,
			     "Updater.Ionization: Must provide elementary charge using 'elemCharge'")
   self._E = assert(tbl.E,
		    "Updater.Ionization: Must provide Voronov constant E using 'E'")
   self._reactRate = tbl.reactRate
   
   if self._reactRate then
      self._A = assert(tbl.A,
		       "Updater.Ionization: Must provide Voronov constant A using 'A'")
      self._K = assert(tbl.K,
		       "Updater.Ionization: Must provide Voronov constant K using 'K'")
      self._P = assert(tbl.P,
		       "Updater.Ionization: Must provide Voronov constant P using 'P'")
      self._X = assert(tbl.X,
		       "Updater.Ionization: Must provide Voronov constant X using 'X'")
   end
      
   -- Dimension of phase space.
   self._pDim = self._phaseBasis:ndim()      
   -- Dimension of configuration space.
   self._cDim = self._confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = self._confBasis:id()
   self._polyOrder = self._confBasis:polyOrder()
   self.idxP       = Lin.IntVec(self._pDim)

   -- Number of basis functions.
   self._numBasisC = self._confBasis:numBasis()

   -- Define ionization temperature calculation
   self._IonizationTempCalc = IonizationDecl.ionizationTemp(self._basisID, self._cDim, self._polyOrder)
   -- Define Voronov reaction rate
   self._VoronovReactRateCalc = IonizationDecl.voronovCoef(self._basisID, self._cDim, self._polyOrder)

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._tmEvalMom = 0.0
end

function Ionization:ionizationTemp(elcVtSq, vtSqIz)
   local tmEvalMomStart = Time.clock()
   local grid           = self._onGrid
   
   local confIndexer = elcVtSq:genIndexer()
   local elcVtSqItr  = elcVtSq:get(1)
   local vtSqIzItr   = vtSqIz:get(1)

   local confRange = elcVtSq:localRange()
   if self.onGhosts then confRange = elcVtSq:localExtRange() end
   
   -- Configuration space loop
   for cIdx in confRange:rowMajorIter() do
      grid:setIndex(cIdx)

      elcVtSq:fill(confIndexer(cIdx), elcVtSqItr)
      vtSqIz:fill(confIndexer(cIdx), vtSqIzItr)

      self._IonizationTempCalc(self._elemCharge, self._elcMass, elcVtSqItr:data(), self._E, vtSqIzItr:data())
     
   end
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function Ionization:reactRateCoef(neutM0, neutVtSq, elcVtSq, coefIz) --, cflRateByCell)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local vDim = self._pDim - self._cDim

   local confIndexer = elcVtSq:genIndexer()
   local neutM0Itr   = neutM0:get(1)
   local neutVtSqItr = neutVtSq:get(1)
   local elcVtSqItr  = elcVtSq:get(1)
   local coefIzItr   = coefIz:get(1)

   local confRange = elcVtSq:localRange()
   if self.onGhosts then confRange = elcVtSq:localExtRange() end

   -- Get the interface for setting global CFL frequencies
   local cflRateByCell     = self._cflRateByCell
   local cflRateByCellIdxr = cflRateByCell:genIndexer()
   local cflRateByCellPtr  = cflRateByCell:get(1)
   local cflRange = cflRateByCell:localRange()
   local velRange = cflRange:selectLast(vDim)
   
   local cflRate
   -- Configuration space loop
   for cIdx in confRange:rowMajorIter() do
      grid:setIndex(cIdx)

      neutM0:fill(confIndexer(cIdx), neutM0Itr)
      neutVtSq:fill(confIndexer(cIdx), neutVtSqItr)
      elcVtSq:fill(confIndexer(cIdx), elcVtSqItr)
      coefIz:fill(confIndexer(cIdx), coefIzItr)

      cflRate = self._VoronovReactRateCalc(self._elemCharge, self._elcMass, neutM0Itr:data(), neutVtSqItr:data(), elcVtSqItr:data(), self._E, self._A, self._K, self._P, self._X, coefIzItr:data())
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

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function Ionization:_advance(tCurr, inFld, outFld)

   if not self._reactRate then
      local elcVtSq = inFld[1]
      local vtSqIz  = outFld[1]
      self:ionizationTemp(elcVtSq, vtSqIz)
   else
      --local cflRateByCell = self._cflRateByCell
      local neutM0   = inFld[1]
      local neutVtSq = inFld[2]
      local elcVtSq  = inFld[3]
      local coefIz   = outFld[1]
      self:reactRateCoef(neutM0, neutVtSq, elcVtSq, coefIz)
   end
   
end

function Ionization:evalMomTime() return self._tmEvalMom end

return Ionization
