-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate Voronov ionization for various species
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"

-- BGK Collisions updater object
local VoronovIonization = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function VoronovIonization:init(tbl)
   VoronovIonization.super.init(self, tbl) -- setup base object

   self._confGrid = assert(tbl.confGrid,
			   "Updater.VoronovIonization: Must provide configuration space grid object using 'confGrid'")
   self._confBasis = assert(tbl.confBasis,
			    "Updater.VoronovIonization: Must provide configuration space basis object using 'confBasis'")
   self._phaseGrid = assert(tbl.phaseGrid,
			    "Updater.VoronovIonization: Must provide phase space grid object using 'phaseGrid'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.VoronovIonization: Must provide phase space basis object using 'phaseBasis'")

   self._elcMass = assert(tbl.elcMass,
			  "Updater.VoronovIonization: Must provide electron mass using 'elcMass'")
   self._elemCharge = assert(tbl.elemCharge,
			    "Updater.VoronovIonization: Must provide elementary charge using 'elemCharge'")
   self._A = assert(tbl.A,
		    "Updater.VoronovIonization: Must provide Voronov constant A using 'A'")
   self._E = assert(tbl.E,
		    "Updater.VoronovIonization: Must provide Voronov constant E using 'E'")
   self._K = assert(tbl.K,
		    "Updater.VoronovIonization: Must provide Voronov constant K using 'K'")
   self._P = assert(tbl.P,
		    "Updater.VoronovIonization: Must provide Voronov constant P using 'P'")
   self._X = assert(tbl.X,
		    "Updater.VoronovIonization: Must provide Voronov constant X using 'X'")


   -- Number of quadrature points in each direction
   self._N = tbl.numConfQuad and tbl.numConfQuad or self._confBasis:polyOrder() + 1

   -- 1D configuration space weights and ordinates
   local ordinates = GaussQuadRules.ordinates[self._N]
   local weights = GaussQuadRules.weights[self._N]

   local numConfDims = self._confBasis:ndim()
   local l, u = {}, {}
   for d = 1, numConfDims do l[d], u[d] = 1, self._N end
   local quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   local numOrdinates = quadRange:volume() -- number of ordinates
   -- construct weights and ordinates for integration in multiple dimensions
   self._confOrdinates = Lin.Mat(numOrdinates, numConfDims)
   local nodeNum = 1
   for idx in quadRange:colMajorIter() do
      for d = 1, numConfDims do
	 self._confOrdinates[nodeNum][d] = ordinates[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   local numBasis = self._confBasis:numBasis()
   self._confBasisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
   -- pre-compute values of basis functions at quadrature nodes
   for n = 1, numOrdinates do
      self._confBasis:evalBasis(self._confOrdinates[n],
				self._confBasisAtOrdinates[n])
   end

   local numPhaseDims = self._phaseBasis:ndim()
   for d = 1, numPhaseDims do l[d], u[d] = 1, self._N end
   quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   numOrdinates = quadRange:volume() -- number of ordinates
   -- construct weights and ordinates for integration in multiple dimensions
   self._phaseOrdinates = Lin.Mat(numOrdinates, numPhaseDims)
   self._phaseWeights = Lin.Vec(numOrdinates)
   nodeNum = 1
   for idx in quadRange:colMajorIter() do
      self._phaseWeights[nodeNum] = 1.0
      for d = 1, numPhaseDims do
	 self._phaseWeights[nodeNum] =
	    self._phaseWeights[nodeNum] * weights[idx[d]]
	 self._phaseOrdinates[nodeNum][d] = ordinates[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   numBasis = self._phaseBasis:numBasis()
   self._phaseBasisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
   -- pre-compute values of basis functions at quadrature nodes
   for n = 1, numOrdinates do
      self._phaseBasis:evalBasis(self._phaseOrdinates[n],
				 self._phaseBasisAtOrdinates[n])
   end

   -- timings
   self._tmEvalMom = 0.0
   self._tmProjectMaxwell = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function VoronovIonization:_advance(tCurr, dt, inFld, outFld)
   local numConfDims = self._confGrid:ndim()
   local numConfBasis = self._confBasis:numBasis()
   local numPhaseDims = self._phaseGrid:ndim()
   local numPhaseBasis = self._phaseBasis:numBasis()
   local numVelDims = numPhaseDims - numConfDims

   -- Define regions for looping over ordinates
   local l, u = {}, {}
   for d = 1, numConfDims do l[d], u[d] = 1, self._N end
   local confQuadRange = Range.Range(l, u)
   local confQuadIndexer = Range.makeColMajorGenIndexer(confQuadRange)
   for d = 1, numPhaseDims do l[d], u[d] = 1, self._N end
   local phaseQuadRange = Range.Range(l, u)
   local phaseQuadIndexer = Range.makeColMajorGenIndexer(phaseQuadRange)
   local phaseIdx = {}

   -- Vorozon coefficient at quadrature points
   local numConfOrdinates = confQuadRange:volume()
   local vorozonCoefOrd = Lin.Vec(numConfOrdinates)
   -- Additional variables which we don't want to allocate inside the
   -- loops
   local elcM0Ord, elcM2Ord = 0.0, 0.0
   local elcM1iOrd = Lin.Vec(numVelDims)
   local U, u2, vth2, Te = 0.0, 0.0, 0.0, 0.0
   local confMu, phaseMu = 0, 0
   local offset = 0
   local numPhaseOrdinates = phaseQuadRange:volume()
   local RHS = Lin.Vec(numPhaseBasis)
   local fOrd = Lin.Vec(numPhaseOrdinates)

   -- Timings
   -- local tmEvalMomStart = 0.0
   -- local tmProjectMaxwellStart = 0.0

   -- Get the inputs and outputs
   local elcM0 = assert(inFld[1],
			"VoronovIonization.advance: Must specify an input fluid moments field")
   local elcM1i = assert(inFld[2],
			 "VoronovIonization.advance: Must specify an input fluid moments field")
   local elcM2 = assert(inFld[3],
			"VoronovIonization.advance: Must specify an input fluid moments field")

   local fElcOut = assert(outFld['elc'],
			  "VoronovIonization.advance: Must specify an electron output field")
   local fIonOut = assert(outFld['ion'],
			  "VoronovIonization.advance: Must specify an ion output field")
   local fNeutOnElcOut = assert(outFld['neutOnElc'],
				"VoronovIonization.advance: Must specify a neutral output field on electron mesh")
   local fNeutOnIonOut = assert(outFld['neutOnIon'],
				"VoronovIonization.advance: Must specify a neutral output field on ion mesh")

   local elcM0Itr = elcM0:get(1)
   local elcM1iItr = elcM1i:get(1)
   local elcM2Itr = elcM2:get(1)
   local fElcItr = fElcOut:get(1)
   local fIonItr = fIonOut:get(1)
   local fNeutOnElcItr = fNeutOnElcOut:get(1)
   local fNeutOnIonItr = fNeutOnIonOut:get(1)

   -- Get the Ranges to loop over the domain
   local confRange = elcM0:localRange()
   local confIndexer = elcM0:genIndexer()
   local phaseRangeElc = fElcOut:localRange()
   local phaseIndexerElc = fElcOut:genIndexer()
   local phaseRangeIon = fIonOut:localRange()
   local phaseIndexerIon = fIonOut:genIndexer()
   local le, ue, li, ui = {}, {}, {}, {}
   for d = 1, numVelDims do
      le[d] = phaseRangeElc:lower(numConfDims + d)
      ue[d] = phaseRangeElc:upper(numConfDims + d)
      li[d] = phaseRangeIon:lower(numConfDims + d)
      ui[d] = phaseRangeIon:upper(numConfDims + d)
   end
   local velRangeElc = Range.Range(le, ue)
   local velRangeIon = Range.Range(li, ui)

   -- Configuration space loop
   for confIdx in confRange:colMajorIter() do
      elcM0:fill(confIndexer(confIdx), elcM0Itr)
      elcM1i:fill(confIndexer(confIdx), elcM1iItr)
      elcM2:fill(confIndexer(confIdx), elcM2Itr)

      -- Evaluate the the moments (given as expansion coefficiens)
      -- on the ordinates
      -- local tmEvalMomStart = Time.clock()
      for muIdx in confQuadRange:colMajorIter() do
	 confMu = confQuadIndexer(muIdx)
	 elcM0Ord = 0
	 for d = 1, numVelDims do elcM1iOrd[d] = 0 end
	 elcM2Ord = 0
	 for k = 1, numConfBasis do
	    elcM0Ord = elcM0Ord +
	       elcM0Itr[k] * self._confBasisAtOrdinates[confMu][k]
	    elcM2Ord = elcM2Ord +
	       elcM2Itr[k] * self._confBasisAtOrdinates[confMu][k]
	 end
	 offset = 0
	 for d = 1, numVelDims do
	    for k = 1, numConfBasis do
	       elcM1iOrd[d] = elcM1iOrd[d] +
		  elcM1iItr[offset + k] * self._confBasisAtOrdinates[confMu][k]
	    end
	    offset = offset + numConfBasis
	 end

	 -- Calculate the temperature
	 u2 = 0
	 for d = 1, numVelDims do
	    u2 = u2 + elcM1iOrd[d]*elcM1iOrd[d] / (elcM0Ord*elcM0Ord)
	 end
	 vth2 = elcM2Ord / elcM0Ord - u2
	 Te = 0.5*self._elcMass*vth2 / self._elemCharge -- temperature in eV
	 U = self._E/Te
	 vorozonCoefOrd[confMu] = elcM0Ord * self._A * (1+self._P*U) / (self._X+U) *
	    math.pow(U, self._K) * math.exp(-U)
      end
      -- self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart

      -- tmProjectMaxwellStart = Time.clock()
      -- Velocity space loop for electrons
      for velIdx in velRangeElc:colMajorIter() do
	 -- Construct the phase space index ot of the configuration
	 -- space a velocity space indices
	 for d = 1, numConfDims do phaseIdx[d] = confIdx[d] end
	 for d = 1, numVelDims do phaseIdx[d + numConfDims] = velIdx[d] end
	 fElcOut:fill(phaseIndexerElc(phaseIdx), fElcItr)
	 fNeutOnElcOut:fill(phaseIndexerElc(phaseIdx), fNeutOnElcItr)

	 for k = 1, numPhaseBasis do RHS[k] = 0 end
	 for muIdx in phaseQuadRange:colMajorIter() do
	    confMu = confQuadIndexer(muIdx)
	    phaseMu = phaseQuadIndexer(muIdx)

	    -- Project on basis
	    fOrd = 0.0
	    for k = 1, numPhaseBasis do -- evaluation
	       fOrd = fOrd + self._phaseBasisAtOrdinates[phaseMu][k] * fNeutOnElcItr[k]
	    end
	    for k = 1, numPhaseBasis do -- itegral
	       RHS[k] = RHS[k] +
		  self._phaseWeights[phaseMu] * vorozonCoefOrd[confMu] *
		  self._phaseBasisAtOrdinates[phaseMu][k] * fOrd
	    end
	 end
	 -- Modify the solutions
	 for k = 1, numPhaseBasis do
	    fElcItr[k] = fElcItr[k] + dt * RHS[k]
	    fNeutOnElcItr[k] = fNeutOnElcItr[k] - dt * RHS[k]
	 end
      end

      -- Velocity space loop for electrons
      for velIdx in velRangeIon:colMajorIter() do
	 -- Construct the phase space index ot of the configuration
	 -- space a velocity space indices
	 for d = 1, numConfDims do phaseIdx[d] = confIdx[d] end
	 for d = 1, numVelDims do phaseIdx[d + numConfDims] = velIdx[d] end
	 fIonOut:fill(phaseIndexerIon(phaseIdx), fIonItr)
	 fNeutOnIonOut:fill(phaseIndexerIon(phaseIdx), fNeutOnIonItr)

	 for k = 1, numPhaseBasis do RHS[k] = 0 end
	 for muIdx in phaseQuadRange:colMajorIter() do
	    confMu = confQuadIndexer(muIdx)
	    phaseMu = phaseQuadIndexer(muIdx)

	    -- Project on basis
	    fOrd = 0.0
	    for k = 1, numPhaseBasis do -- evaluation
	       fOrd = fOrd + self._phaseBasisAtOrdinates[phaseMu][k] * fNeutOnIonItr[k]
	    end
	    for k = 1, numPhaseBasis do -- itegral
	       RHS[k] = RHS[k] +
		  self._phaseWeights[phaseMu] * vorozonCoefOrd[confMu] *
		  self._phaseBasisAtOrdinates[phaseMu][k] * fOrd
	    end
	 end
	 -- Modify the solutions
	 for k = 1, numPhaseBasis do
	    fIonItr[k] = fIonItr[k] + dt * RHS[k]
	    fNeutOnIonItr[k] = fNeutOnIonItr[k] - dt * RHS[k]
	 end
      end
      -- self._tmProjectMaxwell = self._tmProjectMaxwell + Time.clock() -
      -- 	 tmProjectMaxwellStart
   end
   return true, GKYL_MAX_DOUBLE
end

function VoronovIonization:evalMomTime() return self._tmEvalMom end
function VoronovIonization:projectMaxwellTime() return self._tmProjectMaxwell end

return VoronovIonization
