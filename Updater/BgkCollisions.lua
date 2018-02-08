-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate BGK collisions for various species
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"

-- BGK Collisions updater object
local BgkCollisions = Proto(UpdaterBase)

function BgkCollisions:init(tbl)
   BgkCollisions.super.init(self, tbl) -- setup base object

   self._isFirst = true -- for use in advance()

   self._speciesList = assert(tbl.speciesList, "Updater.BgkCollisions: Must provide a list of species using 'speciesList'")
   self._collFreq = assert(tbl.collFreq, "Updater.BgkCollisions: Must provide a collision frequency using 'collFreq'")

   self._confGrid = assert(tbl.confGrid, "Updater.BgkCollisions: Must provide configuration space grid object using 'confGrid'")
   self._confBasis = assert(tbl.confBasis, "Updater.BgkCollisions: Must provide configuration space basis object using 'confBasis'")
   self._phaseGrid = assert(tbl.phaseGrid, "Updater.BgkCollisions: Must provide phase space grid object using 'phaseGrid'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.BgkCollisions: Must provide phase space basis object using 'phaseBasis'")

   assert(#self._speciesList == #self._collFreq, "List of species table and collision frequences table must have the same size")
   assert(#self._speciesList == #self._phaseGrid, "List of species table and phase space grid table must have the same size")
   assert(#self._speciesList == #self._phaseBasis, "List of species table and phase space basis table must have the same size")

   for _, nm in ipairs(speciesList) do
      local N = tbl.numQuad[mn] and tbl.numQuad[mn] or self._phaseBasis[mn]:polyOrder()+1 -- number of quadrature points in each direction

      -- 1D weights and ordinates
      local ordinates = GaussQuadRules.ordinates[N],
      local weights = GaussQuadRules.weights[N]

      local ndim = self._phasePasis[nm]:ndim()
      local l, u = {}, {}
      for d = 1, ndim do l[d], u[d] = 1, N end
      local quadRange = Range.Range(l, u) -- for looping over quadrature nodes
      local numOrdinates = quadRange:volume() -- number of ordinates
   
      -- construct weights and ordinates for integration in multiple dimensions
      self._ordinates = Lin.Mat(numOrdinates, ndim)
      self._weights = Lin.Vec(numOrdinates)
      local nodeNum = 1
      for idx in quadRange:colMajorIter() do
	 self._weights[nodeNum] = 1.0
	 for d = 1, ndim do
	    self._weights[nodeNum] = self._weights[nodeNum]*weights[idx[d]]
	    self._ordinates[nodeNum][d] = ordinates[idx[d]]
	 end
	 nodeNum = nodeNum + 1
      end

      local numBasis = self._basis:numBasis()
      self._basisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
      -- pre-compute values of basis functions at quadrature nodes
      for n = 1, numOrdinates do
	 self._basis:evalBasis(self._ordinates[n], self._basisAtOrdinates[n])
      end
      
      -- construct various functions from template representations
      self._compToPhys = loadstring(compToPhysTempl {NDIM = ndim} )()
   end
end

-- advance method
function BgkCollisions:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local qOut = assert(outFld[1], "BgkCollisions.advance: Must specify an output field")

   local ndim = grid:ndim()
   local numOrd = #self._weights
   local numBasis = self._basis:numBasis()
   local numVal = qOut:numComponents()/numBasis

   if self._isFirst then
      -- construct function to evaluate function at specified coorindate
      self._evalFunc = loadstring(evalFuncTempl { M = numVal } )()
   end

   -- sanity check: ensure number of variables, components and basis functions are consistent
   assert(qOut:numComponents() % numBasis == 0, "BgkCollisions:advance: Incompatible input field")

   local dx = Lin.Vec(ndim) -- cell shape
   local xc = Lin.Vec(ndim) -- cell center
   local fv = Lin.Mat(numOrd, numVal) -- function values at ordinates
   local xmu = Lin.Vec(ndim) -- coordinate at ordinate

   local localRange = qOut:localRange()
   local indexer = qOut:genIndexer() 
   local fItr = qOut:get(1)

   -- loop, computing projections in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      -- get cell shape, cell center coordinates
      for d = 1, ndim do dx[d] = grid:dx(d) end
      grid:cellCenter(xc)

      -- precompute value of function at each ordinate
      for mu = 1, numOrd do
	 self._compToPhys(self._ordinates[mu], dx, xc, xmu) -- compute coordinate
	 self._evalFunc(tCurr, xmu, self._evaluate, fv[mu]) 
      end

      qOut:fill(indexer(idx), fItr)

      local offset = 0
      -- update each component
      for n = 1, numVal do
	 -- update each expansion coefficient
	 for k = 1, numBasis do
	    fItr[offset+k] = 0.0
	    -- loop over quadrature points, accumulating contribution to expansion coefficient
	    for mu = 1, numOrd do
	       fItr[offset+k] = fItr[offset+k] + self._weights[mu]*self._basisAtOrdinates[mu][k]*fv[mu][n]
	    end
	 end
	 offset = offset+numBasis
      end
   end

   self._isFirst = false
   return true, GKYL_MAX_DOUBLE
end

return BgkCollisions
