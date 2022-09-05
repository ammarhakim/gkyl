-- Gkyl ------------------------------------------------------------------------
--
-- Updater for projecting flux of distf onto (ghosts).
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase        = require "Updater.Base"
local Proto              = require "Lib.Proto"
local ProjectFluxModDecl = require "Updater.projectFluxOnGhostsData.ProjectFluxModDecl"
local xsys               = require "xsys"
local Lin                = require "Lib.Linalg"
local Time               = require "Lib.Time"

-- Project flux onto ghosts updater object.
local ProjectFluxFunc = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function ProjectFluxFunc:init(tbl)
   ProjectFluxFunc.super.init(self, tbl) -- setup base object
   
   self._onGrid = assert(tbl.onGrid,
                         "Updater.ProjectFluxFunc: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis,
                            "Updater.ProjectFluxFunc: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
                             "Updater.ProjectFluxFunc: Must provide phase space basis object using 'phaseBasis'")
   self._edgeVal = assert(tbl.edgeValue,
                          "Updater.ProjectFluxFunc: Must provide edge values using 'edgeValue' ") 
   self._dir = assert(tbl.direction,
                      "Updater.ProjectFluxFunc: Must provide edge values using 'direction' ") 
   -- Dimension of phase space.
   self._pDim = self._phaseBasis:ndim()      
   -- Dimension of configuration space.
   self._cDim = self._confBasis:ndim()
   -- Dimension of velocity space.
   self._vDim = self._pDim - self._cDim
   
   -- Basis name and polynomial order.
   self._basisID   = self._confBasis:id()
   self._polyOrder = self._confBasis:polyOrder()
   self.idxP       = Lin.IntVec(self._pDim)

   -- Number of basis functions.
   self._numBasisC = self._confBasis:numBasis()

   -- Define ionization temperature calculation
   self._projectFlux = ProjectFluxModDecl.projectFlux(self._basisID, self._cDim, self._vDim, self._dir, self._polyOrder)

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)
end

function ProjectFluxFunc:_advance(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()
   local grid           = self._onGrid
   local cDim, vDim, pDim = self._cDim, self._vDim, self._pDim

   -- Get vpar limits of cell.
   local vpardir = cDim*2 -- (vx for 1x3v, vz for 3x3v)
   
   local fIn, fHat = inFld[1], outFld[1]
   fHat:clear(0.0)

   local phaseRange = fIn:localRange()
   if self.onGhosts then -- Extend range to config-space ghosts.
      for dir = 1, cDim do 
         phaseRange = phaseRange:extendDir(dir, fIn:lowerGhost(), fIn:upperGhost())
      end
   end

   local fInItr, fHatItr = fIn:get(1), fHat:get(1)

   -- Construct ranges for nested loops.
   local confRange = phaseRange:selectFirst(cDim)
   local velRange  = phaseRange:selectLast(vDim)

   local phaseIndexer = fIn:genIndexer()

   -- Outer loop is threaded and over configuration space.
   for cIdx in confRange:rowMajorIter() do
      
      cIdx:copyInto(self.idxP)
   
      -- Inner loop is over velocity space: no threading to avoid race conditions.
      for vIdx in velRange:rowMajorIter() do

	 for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end
	 grid:setIndex(self.idxP)
	 local wv = grid:cellCenterInDir(vpardir)
	 local dv = grid:dx(vpardir)

	 fIn:fill(phaseIndexer(self.idxP), fInItr)
	 fHat:fill(phaseIndexer(self.idxP), fHatItr)

	 self._projectFlux(wv, dv, self._edgeVal, fInItr:data(), fHatItr:data())
      end
   end
end

return ProjectFluxFunc
