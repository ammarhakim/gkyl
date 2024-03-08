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

function VlasovEmissionBC:Maxwellian(t, xn, vt)
   local vSq = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
      vSq = vSq + xn[d]^2
   end
   return math.exp(-vSq/(2.0*vt^2))
end

function VlasovEmissionBC:FurmanPiviElastic(t, xn, P1_inf, P1_hat, E_hat, W, p)
   local E = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
      E = E + 0.5*self.mass*xn[d]^2/math.abs(self.charge)
   end
   return P1_inf + (P1_hat - P1_inf)*math.exp((-math.abs(E - E_hat)/W)^p/p)
end

function VlasovEmissionBC:CazauxElastic(t, xn, E_f, phi)
   local E = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
      E = E + 0.5*self.mass*xn[d]^2/math.abs(self.charge)
   end
   local E_s = E + E_f + phi 
   local G = 1 + (E_s - E)/E
   return (1 - math.sqrt(G))^2/(1 + math.sqrt(G))^2
end

function VlasovEmissionBC:ConstantElastic(t, xn, gain)
   return gain
end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VlasovEmissionBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.saveFlux = tbl.saveFlux or false

   self.bcKind = assert(tbl.spectrum, "VlasovEmissionBC: must specify the type of emission spectrum in 'bcKind'.")
   self.gammaKind = assert(tbl.yield, "VlasovEmissionBC: must specify the type of emission spectrum in 'bcKind'.")
   self.inSpecies = assert(tbl.inSpecies, "VlasovEmissionBC: must specify names of impacting species in 'inSpecies'.")
   self.elasticKind = tbl.elastic or "none"
   self.bcParam = {}
   self.gammaParam = {}
   self.proj = {}
   self.mass = mySpecies.mass
   self.charge = mySpecies.charge
   self.tbound = tbl.tbound
   for ispec, otherNm in ipairs(self.inSpecies) do
      self.bcParam[otherNm] = Lin.Vec(10)
      self.bcParam[otherNm]:data()[0] = self.mass
      self.bcParam[otherNm]:data()[1] = self.charge
      if self.bcKind == "chung-everhart" then
	 self.bcParam[otherNm]:data()[2] = assert(tbl.spectrumFit[ispec].phi, "VlasovEmissionBC: must specify the material work function in 'phi'.")
         self.proj[otherNm] = Projection.VlasovProjection.FunctionProjection
	    { func = function(t, zn) return self:ChungEverhart(t, zn, self.bcParam[otherNm]:data()[2]) end, }
      elseif self.bcKind == "gaussian" then
         self.bcParam[otherNm]:data()[2] = assert(tbl.spectrumFit[ispec].E0, "VlasovEmissionBC: must specify fitting parameter 'E0'.")
         self.bcParam[otherNm]:data()[3] = assert(tbl.spectrumFit[ispec].tau, "VlasovEmissionBC: must specify fitting parameter 'tau'.")
	 self.proj[otherNm] = Projection.VlasovProjection.FunctionProjection
	    { func = function(t, zn) return self:Gaussian(t, zn, self.bcParam[otherNm]:data()[2], self.bcParam[otherNm]:data()[3]) end, }
      elseif self.bcKind == "maxwellian" then
	 self.bcParam[otherNm]:data()[2] = assert(tbl.spectrumFit[ispec].vt, "VlasovEmissionBC: must specify the material work function in 'vt'.")
         self.proj[otherNm] = Projection.VlasovProjection.FunctionProjection
	    { func = function(t, zn) return self:Maxwellian(t, zn, self.bcParam[otherNm]:data()[2]) end, }
      else
         assert(false, "VlasovEmissionBC: Fitting model not recognized.")
      end
      self.gammaParam[otherNm] = Lin.Vec(10)
      if self.gammaKind == "furman-pivi" then
	 self.gammaParam[otherNm]:data()[0] = assert(tbl.yieldFit[ispec].mass, "VlasovEmissionBC: must specify the impacting species mass in 'mass'.")
	 self.gammaParam[otherNm]:data()[1] = assert(tbl.yieldFit[ispec].charge, "VlasovEmissionBC: must specify the impacting species charge in 'charge'.")
	 self.gammaParam[otherNm]:data()[2] = assert(tbl.yieldFit[ispec].gammahat_ts, "VlasovEmissionBC: must specify fitting parameter 'gammahat_ts'.")
	 self.gammaParam[otherNm]:data()[3] = assert(tbl.yieldFit[ispec].Ehat_ts, "VlasovEmissionBC: must specify fitting parameter 'Ehat_ts'.")
	 self.gammaParam[otherNm]:data()[4] = assert(tbl.yieldFit[ispec].t1, "VlasovEmissionBC: must specify fitting parameter 't1'.")
	 self.gammaParam[otherNm]:data()[5] = assert(tbl.yieldFit[ispec].t2, "VlasovEmissionBC: must specify fitting parameter 't2'.")
	 self.gammaParam[otherNm]:data()[6] = assert(tbl.yieldFit[ispec].t3, "VlasovEmissionBC: must specify fitting parameter 't3'.")
	 self.gammaParam[otherNm]:data()[7] = assert(tbl.yieldFit[ispec].t4, "VlasovEmissionBC: must specify fitting parameter 't4'.")
	 self.gammaParam[otherNm]:data()[8] = assert(tbl.yieldFit[ispec].s, "VlasovEmissionBC: must specify fitting parameter 's'.")
      elseif self.gammaKind == "schou" then
	 self.gammaParam[otherNm]:data()[0] = assert(tbl.yieldFit[ispec].mass, "VlasovEmissionBC: must specify the impacting species mass in 'mass'.")
	 self.gammaParam[otherNm]:data()[1] = assert(tbl.yieldFit[ispec].charge, "VlasovEmissionBC: must specify the impacting species charge in 'charge'.")
	 self.gammaParam[otherNm]:data()[2] = assert(tbl.yieldFit[ispec].intWall, "VlasovEmissionBC: must specify fitting parameter 'intWall'.")
	 self.gammaParam[otherNm]:data()[3] = assert(tbl.yieldFit[ispec].A2, "VlasovEmissionBC: must specify fitting parameter 'A2'.")
	 self.gammaParam[otherNm]:data()[4] = assert(tbl.yieldFit[ispec].A3, "VlasovEmissionBC: must specify fitting parameter 'A3'.")
	 self.gammaParam[otherNm]:data()[5] = assert(tbl.yieldFit[ispec].A4, "VlasovEmissionBC: must specify fitting parameter 'A4'.")
	 self.gammaParam[otherNm]:data()[6] = assert(tbl.yieldFit[ispec].A5, "VlasovEmissionBC: must specify fitting parameter 'A5'.")
	 self.gammaParam[otherNm]:data()[7] = assert(tbl.yieldFit[ispec].nw, "VlasovEmissionBC: must specify fitting parameter 'nw'.")
      elseif self.gammaKind == "constant" then
	 self.gammaParam[otherNm]:data()[0] = assert(tbl.yieldFit[ispec].mass, "VlasovEmissionBC: must specify the impacting species mass in 'mass'.")
	 self.gammaParam[otherNm]:data()[1] = assert(tbl.yieldFit[ispec].charge, "VlasovEmissionBC: must specify the impacting species charge in 'charge'.")
	 self.gammaParam[otherNm]:data()[2] = assert(tbl.yieldFit[ispec].gain, "VlasovEmissionBC: must specify fitting parameter 'gain'.")
      else
         assert(false, "VlasovEmissionBC: SEY model not recognized.")   
      end
   end

   self.elasticParam = Lin.Vec(10)
   if self.elasticKind == "furman-pivi" then
      self.elastic = true
      self.elasticParam:data()[0] = assert(tbl.elasticFit.P1_inf, "VlasovEmissionBC: must specify fitting parameter 'P1_inf'.")
      self.elasticParam:data()[1] = assert(tbl.elasticFit.P1_hat, "VlasovEmissionBC: must specify fitting parameter 'P1_hat'.")
      self.elasticParam:data()[2] = assert(tbl.elasticFit.E_hat, "VlasovEmissionBC: must specify fitting parameter 'E_hat'.")
      self.elasticParam:data()[3] = assert(tbl.elasticFit.W, "VlasovEmissionBC: must specify fitting parameter 'W'.")
      self.elasticParam:data()[4] = assert(tbl.elasticFit.p, "VlasovEmissionBC: must specify fitting parameter 'p'.")
      self.elasticProj = Projection.VlasovProjection.FunctionProjection
         { func = function(t, zn) return self:FurmanPiviElastic(t, zn, self.elasticParam:data()[0], self.elasticParam:data()[1], self.elasticParam:data()[2], self.elasticParam:data()[3], self.elasticParam:data()[4]) end, }
   elseif self.elasticKind == "cazaux" then
      self.elastic = true
      self.elasticParam:data()[0] = assert(tbl.elasticFit.E_f, "VlasovEmissionBC: must specify fitting parameter 'E_f'.")
      self.elasticParam:data()[1] = assert(tbl.elasticFit.phi, "VlasovEmissionBC: must specify fitting parameter 'phi'.")
      self.elasticProj = Projection.VlasovProjection.FunctionProjection
         { func = function(t, zn) return self:CazauxElastic(t, zn, self.elasticParam:data()[0], self.elasticParam:data()[1]) end, }
   elseif self.elasticKind == "sydorenko" then
      self.elastic = true
      self.elasticParam:data()[0] = assert(tbl.elasticFit.delta_m, "VlasovEmissionBC: must specify fitting parameter 'delta_m'.")
      self.elasticParam:data()[1] = assert(tbl.elasticFit.E_m, "VlasovEmissionBC: must specify fitting parameter 'E_m'.")
      self.elasticParam:data()[2] = assert(tbl.elasticFit.E_0, "VlasovEmissionBC: must specify fitting parameter 'E_0'.")
      self.elasticParam:data()[3] = assert(tbl.elasticFit.dE, "VlasovEmissionBC: must specify fitting parameter 'dE'.")
      self.elasticProj = Projection.VlasovProjection.FunctionProjection
         { func = function(t, zn) return self:SydorenkoElastic(t, zn, self.elasticParam:data()[0], self.elasticParam:data()[1], self.elasticParam:data()[2], self.elasticParam:data()[3]) end, }
   elseif self.elasticKind == "constant" then
      self.elastic = true
      self.elasticParam:data()[0] = assert(tbl.elasticFit.gain, "VlasovEmissionBC: must specify fitting parameter 'gain'.")
      self.elasticProj = Projection.VlasovProjection.FunctionProjection
         { func = function(t, zn) return self:ConstantElastic(t, zn, self.elasticParam:data()[0]) end, }
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
   if self.elastic then
      self.elasticProj:fullInit(species[self.speciesName])
      self.elasticProj:createSolver(species[self.speciesName])
   end
   self.fluxBC = {}
   for _, otherNm in ipairs(self.inSpecies) do
      self.proj[otherNm]:fullInit(species[self.speciesName])
      self.proj[otherNm]:createSolver(species[self.speciesName])
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
   self.globalGhostRange = self.bcEdge=="lower" and distf:localGhostRangeLower()[self.bcDir]
                                                  or distf:localGhostRangeUpper()[self.bcDir]
   self:createBoundaryGrid(self.globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())
   -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
   self:createConfBoundaryGrid(self.globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())
   
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

   self.bcBuffer = allocDistf() -- Buffer used by EmissionSpectrumBc updater.

   if self.elastic then
      self.elasticBuffer1 = allocDistf() -- Buffer used by BasicBc updater.
      self.elasticBuffer2 = allocDistf() -- Buffer used by CartFieldBinOp updater.
      self.elasticSolver = Updater.BasicBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   basis  = self.basis,
         dir     = self.bcDir,  bcType = "reflect",
         onField = mySpecies:rkStepperFields()[1],
      }
      self.weakMultiply = Updater.CartFieldBinOp {
	 weakBasis = self.basis, operation = "Multiply",
	 onGhosts = true,
      }
   end
   
   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then

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
	 model_id   = "GKYL_MODEL_DEFAULT",
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
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0}, self.bcBuffer:getMetaData())
   end

   if self.elastic then
      self.fElasticProj = allocDistf()
      self.elasticProj:advance(0.0, {}, {self.fElasticProj})
   end 

   local mySpecies = species[self.speciesName]

   local distf = mySpecies:getDistF()

   self.bcSolver = {}
   self.negRange = {}
   self.posRange = {}
   self.gamma = {}
   for ispec, otherNm in ipairs(self.inSpecies) do
      local bc = self.fluxBC[otherNm]
      self.gamma[otherNm] = self:allocCartField(bc.boundaryGrid, 1, {0,0}, distf:getMetaData())
      self.bcSolver[otherNm] = Updater.EmissionSpectrumBc{
         onGrid  = bc.boundaryGrid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   vdim   = self.vdim,
         dir     = self.bcDir,  bcType = self.bcKind,
	 gammaType = self.gammaKind,
	 bcParam = self.bcParam[otherNm]:data(), gammaParam = self.gammaParam[otherNm]:data(),
	 onField = self.gamma[otherNm],
      }
      self.negRange[otherNm] = self.bcSolver[otherNm].negRange[self.bcDir]
      self.posRange[otherNm] = self.bcSolver[otherNm].posRange[self.bcDir]
      self.integNumDensityCalc = Updater.DistFuncMomentDG {
         onGrid     = self.boundaryGrid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,          moment     = "M0",
	 model_id   = "GKYL_MODEL_DEFAULT",
         isIntegrated = true,
      }
   end
   
   self.fProj = {}
   self.bcFlux = {}
   self.weight = {}
   self.k = {}
   self.otherSkinRange = {}
   for _, otherNm in ipairs(self.inSpecies) do
      local otherSpecies = species[otherNm]
      self.fProj[otherNm] = allocDistf()
      self.proj[otherNm]:advance(0.0, {}, {self.fProj[otherNm]})
      self.bcFlux[otherNm] = self.fluxBC[otherNm]:allocIntThreeMoments()
      self.weight[otherNm] = self:allocCartField(self.confBoundaryGrid, 2, {0,0}, self.bcBuffer:getMetaData())
      self.k[otherNm] = self:allocCartField(self.confBoundaryGrid, 1, {0,0}, self.bcBuffer:getMetaData())
      self.otherSkinRange[otherNm] = self.bcEdge=="lower" and otherSpecies:getDistF():localSkinRangeLower()[self.bcDir]
	                           or otherSpecies:getDistF():localSkinRangeUpper()[self.bcDir]
   end
end

function VlasovEmissionBC:advanceCrossSpeciesCoupling(tCurr, species, inIdx, outIdx)
   self.bcBuffer:clearRange(0.0, self.bcBuffer:localRange())
   for ispec, otherNm in ipairs(self.inSpecies) do
      local otherSpecies = species[otherNm]
      local bc = self.fluxBC[otherNm]

      local qIn = bc:getBoundaryFluxFields()[outIdx]
      local mout = self.bcFlux[otherNm]

      local confRange = mout:localExtRange()
      local phaseRange = self.bcEdge=="lower" and self.negRange[otherNm] or self.posRange[otherNm]
      bc.integNumDensityCalc:advance(tCurr, {qIn, confRange, phaseRange}, {mout})

      local fOther = otherSpecies:rkStepperFields()[inIdx]
      self.bcSolver[otherNm]:advance(tCurr, {fOther, self.bcParam[bc.speciesName], self.fProj[otherNm], self.bcFlux[otherNm], bc.boundaryGrid, self.gamma[otherNm], self.otherSkinRange[otherNm]}, {self.weight[otherNm], self.k[otherNm], self.bcBuffer})
   end
end

function VlasovEmissionBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   local tScale = 1.0
   if self.tbound and tCurr < self.tbound then
      tScale = math.sin(math.pi*tCurr/(2*self.tbound))
      self.bcBuffer:scale(tScale)
   end
   if self.elastic then
      self.elasticSolver:advance(tCurr, {self.elasticBuffer1}, {fIn})
      self.weakMultiply:advance(tm, {self.elasticBuffer1, self.fElasticProj}, {self.elasticBuffer2})
      self.bcBuffer:accumulateRange(tScale, self.elasticBuffer2, self.elasticBuffer2:localRange())
   end
   fIn:copyRangeToRange(self.bcBuffer, self.globalGhostRange, self.bcBuffer:localRange())
end

function VlasovEmissionBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return VlasovEmissionBC
