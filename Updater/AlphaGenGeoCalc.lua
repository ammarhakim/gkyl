-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate SigmaCX for the Pauls CX model
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local AlphaGenGeoDecl = require "Updater.alphaGenGeoCalcData.AlphaGenGeoModDecl"
local xsys = require "xsys"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

-- Charge exchange collisions updater object.
local AlphaGenGeoCalc = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function AlphaGenGeoCalc:init(tbl)
   AlphaGenGeoCalc.super.init(self, tbl) -- setup base object

   self._confGrid = assert(tbl.confGrid,
			     "Updater.AlphaGenGeo: Must provide conf grid object using 'confGrid'")   
   self._phaseGrid = assert(tbl.onGrid,
			     "Updater.AlphaGenGeo: Must provide phase grid object using 'phaseGrid'")
   self._confBasis = assert(tbl.confBasis,
			     "Updater.AlphaGenGeo: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.AlphaGenGeo: Must provide phase space basis object using 'phaseBasis'")
   
   -- Dimension of spaces.
   self._cDim = self._confBasis:ndim()
   self._pDim = self._phaseBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Basis name and polynomial order.
   self._basisID = self._phaseBasis:id()
   self._polyOrder = self._phaseBasis:polyOrder()

   -- Cell index, center, and dx.
   self.idx = Lin.IntVec(self._pDim)
   self.xc = Lin.Vec(self._pDim)
   self.dx = Lin.Vec(self._pDim)

   -- Number of basis functions.
   self._numBasis = self._confBasis:numBasis()

   self._calcAlphaGenGeo = AlphaGenGeoDecl.AlphaGenGeo(self._basisID, self._cDim, self._vDim, self._polyOrder)
   
   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function AlphaGenGeoCalc:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local confGrid = self._confGrid
   local phaseGrid = self._phaseGrid
   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local tvComp = assert(inFld[1], "AlphaGenGo.advance: Must specify tangent vector components as input[1]")
   local gxx = assert(inFld[2], "AlphaGenGeo.advance: Must specify gij as input[2]-input[7]")
   local gxy = assert(inFld[3], "AlphaGenGeo.advance: Must specify gij as input[2]-input[7]")
   local gxz = assert(inFld[4], "AlphaGenGeo.advance: Must specify gij as input[2]-input[7]")
   local gyy = assert(inFld[5], "AlphaGenGeo.advance: Must specify gij as input[2]-input[7]")
   local gyz = assert(inFld[6], "AlphaGenGeo.advance: Must specify gij as input[2]-input[7]")
   local gzz = assert(inFld[7], "AlphaGenGeo.advance: Must specify gij as input[2]-input[7]")

   local alphaGeo = assert(outFld[1], "AlphaGenGeo.advance: Must specify alphaGeo as output")

   local tvCompItr = tvComp:get(1)
   local gxxItr, gxyItr, gxzItr = gxx:get(1), gxy:get(1), gxz:get(1)
   local gyyItr, gyzItr, gzzItr = gyy:get(1), gyz:get(1), gzz:get(1)

   local alphaItr = alphaGeo:get(1)
   
   local confIndexer = gxx:genIndexer()
   local phaseIndexer = alphaGeo:genIndexer()
   local tvPhaseIndexer = tvComp:genIndexer()
   local phaseRange
   if self.onGhosts then
      phaseRange = alphaGeo:localExtRange()
   else
      phaseRange = alphaGeo:localRange()
   end
   
   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = phaseGrid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId = phaseGrid:subGridSharedId()   -- Local thread ID.
   
   -- Phase space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      confGrid:setIndex(cIdx)

      gxx:fill(confIndexer(cIdx), gxxItr)
      gxy:fill(confIndexer(cIdx), gxyItr)
      gxz:fill(confIndexer(cIdx), gxzItr)
      gyy:fill(confIndexer(cIdx), gyyItr)
      gyz:fill(confIndexer(cIdx), gyzItr)
      gzz:fill(confIndexer(cIdx), gzzItr)
      
      -- The velocity space loop
      for vIdx in velRange:rowMajorIter() do
	 cIdx:copyInto(self.idx)
	 for d = 1, vDim do self.idx[cDim+d] = vIdx[d] end
	 
	 -- Get cell shape, cell center coordinates
         phaseGrid:setIndex(self.idx)
	 phaseGrid:getDx(self.dx)
	 phaseGrid:cellCenter(self.xc)

	 tvComp:fill(tvPhaseIndexer(self.idx), tvCompItr)
	 alphaGeo:fill(phaseIndexer(self.idx), alphaItr)
	 
	 self._calcAlphaGenGeo(self.xc:data(), self.dx:data(), tvCompItr:data(), gxxItr:data(), gxyItr:data(), gxzItr:data(), gyyItr:data(), gyzItr:data(), gzzItr:data(), alphaItr:data())
      end     
   end
   
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function AlphaGenGeoCalc:evalMomTime() return self._tmEvalMom end

return AlphaGenGeoCalc
