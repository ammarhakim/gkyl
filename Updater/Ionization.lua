-- Gkyl ------------------------------------------------------------------------
--
-- Updater for Voronov ionization calculations:
--     temperature (vtSqIz) for electrons or Voronov reaction rate (coefIz)
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local IonizationDecl = require "Updater.ionizationCalcData.IonizationModDecl"
local xsys = require "xsys"
local Lin = require "Lib.Linalg"
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
      
   -- Dimension of configuration space.
   self._cDim = self._confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID = self._confBasis:id()
   self._polyOrder = self._confBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisC = self._confBasis:numBasis()

   -- Define ionization temperature calculation
   self._IonizationTempCalc = IonizationDecl.ionizationTemp(self._basisID, self._cDim, self._polyOrder)
   -- Define Voronov reaction rate
   self._VoronovReactRateCalc = IonizationDecl.voronovCoef(self._basisID, self._cDim, self._polyOrder)

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   self._tmEvalMom = 0.0
end

function Ionization:ionizationTemp(elcVtSq, vtSqIz)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid

   local confIndexer = elcVtSq:genIndexer()
   local elcVtSqItr = elcVtSq:get(1)
   local vtSqIzItr = vtSqIz:get(1)

   local confRange = elcVtSq:localRange()
   if self.onGhosts then confRange = elcVtSq:localExtRange() end

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- Local thread ID.
   
   -- Configuration space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(cIdx)

      elcVtSq:fill(confIndexer(cIdx), elcVtSqItr)
      vtSqIz:fill(confIndexer(cIdx), vtSqIzItr)

      self._IonizationTempCalc(self._elemCharge, self._elcMass, elcVtSqItr:data(), self._E, vtSqIzItr:data())
     
   end
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function Ionization:reactRateCoef(elcVtSq, coefIz)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid

   local confIndexer = elcVtSq:genIndexer()
   local elcVtSqItr = elcVtSq:get(1)
   local coefIzItr = coefIz:get(1)

   local confRange = elcVtSq:localRange()
   if self.onGhosts then confRange = elcVtSq:localExtRange() end

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- Local thread ID.
   
   -- Configuration space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(cIdx)

      elcVtSq:fill(confIndexer(cIdx), elcVtSqItr)
      coefIz:fill(confIndexer(cIdx), coefIzItr)

      self._VoronovReactRateCalc(self._elemCharge, self._elcMass, elcVtSqItr:data(), self._E, self._A, self._K, self._P, self._X, coefIzItr:data())
     
   end
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function Ionization:_advance(tCurr, inFld, outFld)

   local elcVtSq = inFld[1]
   if not self._reactRate then
      local vtSqIz  = outFld[1]
      self:ionizationTemp(elcVtSq, vtSqIz)
   else
      local coefIz = outFld[1]
      self:reactRateCoef(elcVtSq, coefIz)
   end
   
end

function Ionization:evalMomTime() return self._tmEvalMom end

return Ionization
