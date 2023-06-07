-- Gkyl ------------------------------------------------------------------------
--
-- Emitting boundary condition for a Vlasov species.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase     = require "App.BCs.BCsBase"
local DataStruct  = require "DataStruct"
local Updater     = require "Updater"
local Mpi         = require "Comm.Mpi"
local Projection  = require "App.Projection"
local Proto       = require "Lib.Proto"
local Time        = require "Lib.Time"
local Range       = require "Lib.Range"
local Lin         = require "Lib.Linalg"
local Grid        = require "Grid"
local DiagsApp    = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags = require "App.Diagnostics.VlasovDiagnostics"
local xsys        = require "xsys"

local VlasovEmissionBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function VlasovEmissionBC:init(tbl) self.tbl = tbl end

function VlasovEmissionBC:ChungEverhart(t, xn, phi)
   local E = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
      E = E + 0.5*self.mass*xn[d]^2/math.abs(self.charge)
   end
   return E/(E + phi)^4
end

function VlasovEmissionBC:Gaussian(t, xn, E_0, tau)
   local E = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
      E = E + 0.5*self.mass*xn[d]^2/math.abs(self.charge)
   end
   return math.exp(-math.log(E/E_0)^2/(2.0*tau^2))
end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VlasovEmissionBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.saveFlux = tbl.saveFlux or false

   self.bcKind = assert(tbl.kind, "VlasovEmissionBC: must specify the type of BC in 'kind'.")

   if self.bcKind == "gain" then
      self.bcParam = Lin.Vec(1)
      self.bcParam:data()[0] = assert(tbl.gamma, "VlasovEmissionBC: must specify the emission flux ratio in 'gamma'.")
   elseif self.bcKind == "chung" then
      self.inSpecies = assert(tbl.inSpecies, "VlasovEmissionBC: must specify names of impacting species in 'inSpecies'.")
      self.bcParam = {}
      self.gain = {}
      self.elastic = {}
      self.proj = {}
      self.mass = mySpecies.mass
      self.charge = mySpecies.charge
      for ispec, otherNm in ipairs(self.inSpecies) do
	 self.bcParam[otherNm] = Lin.Vec(1)
         self.bcParam[otherNm]:data()[0] = assert(tbl.work[ispec], "VlasovEmissionBC: must specify the material work function in 'work'.")
         self.gain[otherNm] = assert(tbl.gain[ispec], "VlasovEmissionBC: must give gain array in 'gain'")
         self.elastic[otherNm] = tbl.elastic[ispec] or nil
	 self.proj[otherNm] = Projection.KineticProjection.FunctionProjection
	    { func = function(t, zn) return self:ChungEverhart(t, zn, self.bcParam[otherNm]:data()[0]) end, }
      end
   elseif self.bcKind == "gauss" then
      self.inSpecies = assert(tbl.inSpecies, "VlasovEmissionBC: must specify names of impacting species in 'inSpecies'.")
      self.bcParam = {}
      self.gain = {}
      self.elastic = {}
      self.proj = {}
      self.mass = mySpecies.mass
      self.charge = mySpecies.charge
      for ispec, otherNm in ipairs(self.inSpecies) do
	 self.bcParam[otherNm] = Lin.Vec(4)
         self.bcParam[otherNm]:data()[0] = assert(tbl.E0[ispec], "VlasovEmissionBC: must specify fitting parameter 'E0'.")
         self.bcParam[otherNm]:data()[1] = assert(tbl.tau[ispec], "VlasovEmissionBC: must specify fitting parameter 'tau'.")
	 self.bcParam[otherNm]:data()[2] = self.mass
	 self.bcParam[otherNm]:data()[3] = self.charge
         self.gain[otherNm] = assert(tbl.gain[ispec], "VlasovEmissionBC: must give gain array in 'gain'")
         self.elastic[otherNm] = tbl.elastic[ispec] or nil
	 self.proj[otherNm] = Projection.KineticProjection.FunctionProjection
	    { func = function(t, zn) return self:Gaussian(t, zn, self.bcParam[otherNm]:data()[0], self.bcParam[otherNm]:data()[1]) end, }
      end
   end

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function VlasovEmissionBC:setName(nm) self.name = self.speciesName.."_"..nm end

function VlasovEmissionBC:initCrossSpeciesCoupling(species)
   self.fluxBC = {}
   for _, otherNm in ipairs(self.inSpecies) do
      self.proj[otherNm]:fullInit(species[self.speciesName])
      local otherSpecies = species[otherNm]
      self.fluxBC[otherNm] = otherSpecies.nonPeriodicBCs[string.gsub(self.name,self.speciesName.."_","")]
      self.fluxBC[otherNm]:setSaveFlux(true)
   end
end

function VlasovEmissionBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   local globalGhostRange = self.bcEdge=="lower" and distf:localGhostRangeLower()[self.bcDir]
                                                  or distf:localGhostRangeUpper()[self.bcDir]
   self:createBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())
   
   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0}, distf:getMetaData())
   end

   self.bcBuffer = allocDistf() -- Buffer used by BasicBc updater.

   local bcFunc, skinType
   if self.bcKind == "gain" then
      self.bcSolver = Updater.EmissionBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   basis  = self.basis,
         dir     = self.bcDir,  bcType = self.bcKind,
	 bcField = nil,         bcParam = self.bcParam,
         onField = mySpecies:rkStepperFields()[1],
      }
   elseif self.bcKind == "chung" or self.bcKind == "gauss" then
      self.bcSolver = Updater.EmissionSpectrumBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   vdim   = self.vdim,
	 basis  = self.basis,   confBasis  = self.confBasis, 
         dir     = self.bcDir,  bcType = self.bcKind,
	 gain = self.gain,      bcParam = self.bcParam,
	 elastic = self.elastic,
	 onField = mySpecies:rkStepperFields()[1],
      }
   else
      assert(false, "VlasovEmissionBC: BC kind not recognized.")
   end

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then

      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())

      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
      end
      self.allocIntThreeMoments = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.vdim+2, {0,0}, numDensity:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
      end
      self.allocIntMoment = function(self, comp)
         local metaData = {charge = self.charge,  mass = self.mass,}
         local ncomp = comp or 1
         local f = DataStruct.DynVector{numComponents = ncomp,     writeRank = self.confBoundaryGrid:commSet().writeRank,
                                        metaData      = metaData,  comm      = self.confBoundaryGrid:commSet().comm,}
         return f
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.distfInIdxr        = distf:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = allocDistf()
      end
      self.boundaryFluxRate      = allocDistf()
      self.boundaryFluxFieldPrev = allocDistf()

      -- Part of global ghost range this rank owns.
      self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                      or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField = self:allocMoment()
      -- Range spanning ghost cells.
      self.myGlobalConfGhostRange = self.bcEdge=="lower" and numDensity:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                          or numDensity:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      local jacobGeo = externalField.geo and externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo and externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
      end

      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         self.boundaryFluxFields[rkIdx]:copyRangeToRange(qOut, self.boundaryFluxFields[rkIdx]:localRange(), self.myGlobalGhostRange)
      end
      self.copyBoundaryFluxFieldFunc = function(inIdx, outIdx)
         self.boundaryFluxFields[outIdx]:copy(self.boundaryFluxFields[inIdx])
      end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...)
         local args  = {...} -- Package up rest of args as table.
         local nFlds = #args/2
         self.boundaryFluxFields[outIdx]:combine(a, self.boundaryFluxFields[aIdx])
         for i = 1, nFlds do -- Accumulate rest of the fields.
            self.boundaryFluxFields[outIdx]:accumulate(args[2*i-1], self.boundaryFluxFields[args[2*i]])
         end
      end

      -- Number density calculator. Needed regardless of diagnostics (for recycling BCs).
      self.numDensityCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
         phaseBasis = self.basis,         moment     = "M0",
      }

      -- Integrated number density calculator. Needed regardless of diagnostics (for steady state sources).
      self.integNumDensityCalc = Updater.DistFuncMomentDG {
         onGrid     = self.boundaryGrid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,          moment     = "M0",
         isIntegrated = true,
      }

      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         self.calcBoundaryFluxRateFunc = function(dtIn)
            -- Compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt.
            self.boundaryFluxRate:combine( 1.0/dtIn, self.boundaryFluxFields[1],
                                          -1.0/dtIn, self.boundaryFluxFieldPrev)
            self.boundaryFluxFieldPrev:copy(self.boundaryFluxFields[1])
         end
         -- Set up weak multiplication and division operators (for diagnostics).
         self.confWeakMultiply = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "Multiply",
            onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "Divide",
            onRange   = self.confBoundaryField:localRange(),  onGhosts = false,
         }
         self.confWeakDotProduct = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "DotProduct",
            onGhosts  = true,
         }
         -- Volume integral operator (for diagnostics).
         self.volIntegral = {
            scalar = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = 1,
               basis  = self.confBasis,         quantity      = "V",
            },
            vector = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = self.vdim,
               basis  = self.confBasis,         quantity      = "V",
            },
         }
         -- Moment calculators (for diagnostics).
         self.momDensityCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M1i",
         }
         self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M2",
         }
         self.M2ijCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         moment    = "M2ij",
         }
         self.M3iCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         moment    = "M3i",
         }
         self.divideByJacobGeo = self.jacobGeoInv
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeoInv}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
         self.multiplyByJacobGeo = self.jacobGeo
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeo}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function VlasovEmissionBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function VlasovEmissionBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function VlasovEmissionBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function VlasovEmissionBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function VlasovEmissionBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = VlasovDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

-- These are needed to recycle the VlasovDiagnostics with VlasovEmissionBC.
function VlasovEmissionBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                 self.boundaryFluxRate, self.boundaryFluxRate} end
function VlasovEmissionBC:getFlucF() return self.boundaryFluxRate end

function VlasovEmissionBC:createCouplingSolver(species, field, extField)
   local mySpecies = species[self.speciesName]
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0}, self.bcBuffer:getMetaData())
   end
   
   self.fProj = {}
   self.bcFlux = {}
   for _, otherNm in ipairs(self.inSpecies) do
      self.fProj[otherNm] = mySpecies:allocDistf()
      self.proj[otherNm]:advance(0.0, {}, {self.fProj[otherNm]})
      self.bcFlux[otherNm] = self.fluxBC[otherNm]:allocIntThreeMoments()
   end
   
   self.localEdgeFlux = 0.0
end

function VlasovEmissionBC:advanceCrossSpeciesCoupling(tCurr, species, inIdx, outIdx)
   local mySpecies = species[self.speciesName]
   local fIn = mySpecies:rkStepperFields()[inIdx]
   self.k = {}

   for ispec, otherNm in ipairs(self.inSpecies) do
      local otherSpecies = species[otherNm]
      local bc = self.fluxBC[otherNm]
      
      bc.integNumDensityCalc:advance(tCurr, {bc:getBoundaryFluxFields()[outIdx]}, {self.bcFlux[otherNm]})

      local flux = self.bcFlux[otherNm]
      self.localEdgeFlux = 0.0

      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(0)
      for idx in flux:localExtRangeIter() do
         flux:fill(fluxIndexer(idx), fluxItr)
         self.localEdgeFlux = self.localEdgeFlux + fluxItr[1]
      end

      local fOther = otherSpecies:rkStepperFields()[inIdx]
      self.k[otherNm] = self.bcSolver:advance(tCurr, {fOther, self.fProj[otherNm], self.bcParam[bc.speciesName], self.localEdgeFlux, bc.boundaryGrid, self.gain[otherNm]}, {fIn})
   end
end

function VlasovEmissionBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   fIn:clearRange(0.0, fIn:localGhostRangeUpper()[self.bcDir])
   for ispec, otherNm in ipairs(self.inSpecies) do
      fIn:accumulateRange(self.k[otherNm], self.fProj[otherNm], fIn:localGhostRangeUpper()[self.bcDir])
   end
end

function VlasovEmissionBC:getBoundaryFluxFields() return self.boundaryFluxFields end

-- ................... Classes meant as aliases to simplify input files ...................... --
local VlasovConstantGainBC = Proto(VlasovEmissionBC)
function VlasovConstantGainBC:fullInit(mySpecies)
   self.tbl.kind  = "gain"
   VlasovConstantGainBC.super.fullInit(self, mySpecies)
end

local VlasovChungEverhartBC = Proto(VlasovEmissionBC)
function VlasovChungEverhartBC:fullInit(mySpecies)
   self.tbl.kind  = "chung"
   VlasovChungEverhartBC.super.fullInit(self, mySpecies)
end

local VlasovGaussianEmissionBC = Proto(VlasovEmissionBC)
function VlasovGaussianEmissionBC:fullInit(mySpecies)
   self.tbl.kind  = "gauss"
   VlasovGaussianEmissionBC.super.fullInit(self, mySpecies)
end
-- ................... End of VlasovEmissionBC alias classes .................... --

return {VlasovConstantGain     = VlasovConstantGainBC,
        VlasovChungEverhart    = VlasovChungEverhartBC,
        VlasovGaussianEmission = VlasovGaussianEmissionBC}
