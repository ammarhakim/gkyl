-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate Voronov ionization reaction rate coef for various species
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase    = require "Updater.Base"
local LinearDecomp   = require "Lib.LinearDecomp"
local Proto          = require "Lib.Proto"
local IonizationDecl = require "Updater.ionizationCalcData.IonizationModDecl"
local xsys           = require "xsys"
local Lin            = require "Lib.Linalg"
local Time           = require "Lib.Time"

-- Voronov Collisions updater object.
local VoronovReactRateCoef = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function VoronovReactRateCoef:init(tbl)
   VoronovReactRateCoef.super.init(self, tbl) -- setup base object

   self._onGrid     = assert(tbl.onGrid,
			  "Updater.VoronovReactRateCoef: Must provide grid object using 'onGrid'")
   self._confBasis  = assert(tbl.confBasis,
			    "Updater.VoronovReactRateCoef: Must provide configuration space basis object using 'confBasis'")
   self._elcMass    = assert(tbl.elcMass,
			  "Updater.VoronovReactRateCoef: Must provide electron mass using 'elcMass'")
   self._elemCharge = assert(tbl.elemCharge,
			     "Updater.VoronovReactRateCoef: Must provide elementary charge using 'elemCharge'")
   self._A = assert(tbl.A,
		    "Updater.VoronovReactRateCoef: Must provide Voronov constant A using 'A'")
   self._E = assert(tbl.E,
		    "Updater.VoronovReactRateCoef: Must provide Voronov constant E using 'E'")
   self._K = assert(tbl.K,
		    "Updater.VoronovReactRateCoef: Must provide Voronov constant K using 'K'")
   self._P = assert(tbl.P,
		    "Updater.VoronovReactRateCoef: Must provide Voronov constant P using 'P'")
   self._X = assert(tbl.X,
		    "Updater.VoronovReactRateCoef: Must provide Voronov constant X using 'X'")

   -- Dimension of configuration space.
   self._cDim = self._confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = self._confBasis:id()
   self._polyOrder = self._confBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisC = self._confBasis:numBasis()

   -- Define Voronov reaction rate
   self._VoronovReactRateCalc = IonizationDecl.voronovCoef(self._basisID, self._cDim, self._polyOrder)

   self.onGhosts = xsys.pickBool(false, tbl.onGhosts)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function VoronovReactRateCoef:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid

   local elcVtSq  = inFld[1]

   local confIndexer = elcVtSq:genIndexer()

   local elcVtSqItr = elcVtSq:get(1)

   local coefIz       = outFld[1]
   local coefIzItr    = coefIz:get(1)

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

function VoronovReactRateCoef:evalMomTime() return self._tmEvalMom end

return VoronovReactRateCoef
