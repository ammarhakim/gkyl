-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate a cell-averaged Maxwellian
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase           = require "Updater.Base"
local Proto                 = require "Lib.Proto"
local CellAveMaxwellianDecl = require "Updater.cellAveMaxwellianCalcData.MaxwellianCellAvModDecl"
local DataStruct            = require "DataStruct"
local xsys                  = require "xsys"
local Lin                   = require "Lib.Linalg"
local Time                  = require "Lib.Time"

-- Charge exchange collisions updater object.
local CellAveMaxwellian = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function CellAveMaxwellian:init(tbl)
   CellAveMaxwellian.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid,
			 "Updater.CellAvMax: Must provide grid object using 'onGrid'")
   self._confGrid = assert(tbl.confGrid,
			   "Updater.CellAvMax: Must provide confGrid object using 'confGrid'")
   self._confBasis = assert(tbl.confBasis,
			     "Updater.CellAvMax: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.CellAvMaxw: Must provide phase space basis object using 'phaseBasis'")
   self._kineticSpecies = assert(tbl.kineticSpecies,
				 "Updater.CellAvMaxwellian: Must provide solver type (Vm or Gk) using 'kineticSpecies'")
   if self._kineticSpecies == 'Gk' then
      self.mass = assert(tbl.gkfacs[1],
			 "Updater.CellAvMaxwellian: Must provide mass object in 'gkfacs'")
      self.bmag = assert(tbl.gkfacs[2],
			 "Updater.CellAvMaxwellian: Must provide bmag object in 'gkfacs'")
      self.bmagItr = self.bmag:get(1)
   end
   
   -- Dimension of spaces.
   self._cDim = self._confBasis:ndim()
   self._pDim = self._phaseBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Basis name and polynomial order.
   self._basisID = self._phaseBasis:id()
   self._polyOrder = self._phaseBasis:polyOrder()

   -- Number of basis functions.
   self._numBasis = self._confBasis:numBasis()

   -- Define cell center
   self.xc  = Lin.Vec(self._pDim)

   -- Define functions
   if self._kineticSpecies == 'Vm' then
      self._calcMax = CellAveMaxwellianDecl.CellAvMax(self._basisID, self._cDim, self._vDim, self._polyOrder)
   else
      self._calcGkMax = CellAveMaxwellianDecl.GkCellAvMax(self._basisID, self._cDim, self._vDim, self._polyOrder)
   end
   
   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function CellAveMaxwellian:vlasov(m0, u, vtSq, fMax)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local pDim = self._pDim
   local vDim = self._vDim
   
   local confIndexer = m0:genIndexer()
   local phaseIndexer = fMax:genIndexer()
   local numConfBasis = self._confBasis:numBasis()

   local uDim = u:numComponents()/numConfBasis
   
   local m0Itr = m0:get(1)
   if vDim == uDim then
      self.uItr = u:get(1)
   elseif uDim < vDim and uDim == 1 then -- GK uPar has been passed
      -- create vector moment for u 
      self.u = DataStruct.Field {
	 onGrid        = self._confGrid,
	 numComponents = self._confBasis:numBasis()*vDim,
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self._phaseBasis:polyOrder(),
	    basisType = self._phaseBasis:id()
	 },
      }
      self.uParItr = u:get(1)
      self.uItr = self.u:get(1)
   end
   local vtSqItr = vtSq:get(1)
   local fMaxItr = fMax:get(1)
   
   local phaseRange = self.onGhosts and fMax:localExtRange() or fMax:localRange()
   
   -- Phase space loop
   for pIdx in phaseRange:rowMajorIter() do
      grid:setIndex(pIdx)
      grid:cellCenter(self.xc)

      m0:fill(confIndexer(pIdx), m0Itr)
      vtSq:fill(confIndexer(pIdx), vtSqItr)
      fMax:fill(phaseIndexer(pIdx), fMaxItr)

      if vDim == uDim then
	 u:fill(confIndexer(pIdx), self.uItr)
      elseif uDim < vDim and uDim == 1 then
	 u:fill(confIndexer(pIdx), self.uParItr)
	 for d = 1, vDim-1 do
	    for k = 1, numConfBasis do
	       self.uItr[numConfBasis*(d-1)+k] = 0.0
	    end
	 end
	 for k = 1, numConfBasis do
	    self.uItr[numConfBasis*(vDim-1)+k] = self.uParItr[k]
	 end
      end
	 
      self._calcMax(self.xc:data(), m0Itr:data(), self.uItr:data(), vtSqItr:data(), fMaxItr:data())     
   end
   
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function CellAveMaxwellian:gyrokinetic(m0, u, vtSq, fMax)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local pDim = self._pDim
   
   local confIndexer = m0:genIndexer()
   local phaseIndexer = fMax:genIndexer()
   local numConfBasis = self._confBasis:numBasis()

   local uDim = u:numComponents()/numConfBasis
   
   local m0Itr = m0:get(1)
   local vtSqItr = vtSq:get(1)
   local fMaxItr = fMax:get(1)

   if uDim == 1 then
      self.uItr = u:get(1)
   else
      self.uItr = m0:get(1)
      self.uInItr = u:get(1)
   end
   
   local phaseRange = self.onGhosts and fMax:localExtRange() or fMax:localRange()
   
   -- Phase space loop
   for pIdx in phaseRange:rowMajorIter() do
      grid:setIndex(pIdx)
      grid:cellCenter(self.xc)

      m0:fill(confIndexer(pIdx), m0Itr)
      vtSq:fill(confIndexer(pIdx), vtSqItr)
      self.bmag:fill(confIndexer(pIdx), self.bmagItr)
      fMax:fill(phaseIndexer(pIdx), fMaxItr)

      if uDim == 1 then
	 u:fill(confIndexer(pIdx), self.uItr)
      else
	 u:fill(confIndexer(pIdx), self.uInItr)
	 for k = 1, numConfBasis do
	    self.uItr[k] = self.uInItr[numConfBasis*(uDim-1)+k]
	 end
      end
      
      self._calcGkMax(self.mass, self.xc:data(), m0Itr:data(), self.uItr:data(), vtSqItr:data(), self.bmagItr:data(), fMaxItr:data())     
   end
   
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function CellAveMaxwellian:_advance(tCurr, inFld, outFld)

   local m0 = assert(inFld[1], "cellAveMaxwellian: Must specify particle density as input[1]")
   local u = assert(inFld[2], "cellAveMaxwellian: Must specify fluid velocity as input[2]")
   local vtSq = assert(inFld[3], "cellAveMaxwellian: Must specify squared thermal velocity as input[3]")
   local fMax = assert(outFld[1], "cellAveMaxwellian: Must specify an output field")
   if self._kineticSpecies == 'Vm' then
      self:vlasov(m0, u, vtSq, fMax)
   elseif self._kineticSpecies == 'Gk' then
      self:gyrokinetic(m0, u, vtSq, fMax)
   end
   
end
   
function CellAveMaxwellian:evalMomTime() return self._tmEvalMom end

return CellAveMaxwellian
