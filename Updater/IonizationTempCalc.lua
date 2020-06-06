-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate ionizaton temperature (vtSqIz) for electrons
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
local IonizationTemp = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function IonizationTemp:init(tbl)
   IonizationTemp.super.init(self, tbl) -- setup base object

   self._onGrid     = assert(tbl.onGrid,
			  "Updater.IonizationTemp: Must provide grid object using 'onGrid'")
   self._confBasis  = assert(tbl.confBasis,
			    "Updater.IonizationTemp: Must provide configuration space basis object using 'confBasis'")
   self._elcMass    = assert(tbl.elcMass,
			  "Updater.IonizationTemp: Must provide electron mass using 'elcMass'")
   self._elemCharge = assert(tbl.elemCharge,
			     "Updater.IonizationTemp: Must provide elementary charge using 'elemCharge'")
   self._E = assert(tbl.E,
		    "Updater.IonizationTemp: Must provide Voronov constant E using 'E'")

   -- Dimension of configuration space.
   self._cDim = self._confBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = self._confBasis:id()
   self._polyOrder = self._confBasis:polyOrder()

   -- Number of basis functions.
   self._numBasisC = self._confBasis:numBasis()

   -- Define ionization temperature calculation
   self._IonizationTempCalc = IonizationDecl.ionizationTemp(self._basisID, self._cDim, self._polyOrder)

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function IonizationTemp:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid

   local elcVtSq  = inFld[1]

   local confIndexer = elcVtSq:genIndexer()

   local elcVtSqItr = elcVtSq:get(1)

   local vtSqIz       = outFld[1]
   local vtSqIzItr    = vtSqIz:get(1)

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

function IonizationTemp:evalMomTime() return self._tmEvalMom end

return IonizationTemp
