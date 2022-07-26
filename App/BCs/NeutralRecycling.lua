-- Gkyl ------------------------------------------------------------------------
--
-- Electron reflection BC based on the Bronold & Fehske model.
--   F.X. Bronold, H. Fehske, Absorption of an electron by a dielectric wall, Phys. Rev. Lett. 115 (22) (2015) 225001.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase        = require "App.BCs.BCsBase"
local DataStruct     = require "DataStruct"
local Updater        = require "Updater"
local Mpi            = require "Comm.Mpi"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Range          = require "Lib.Range"
local Lin            = require "Lib.Linalg"
local LinearDecomp   = require "Lib.LinearDecomp"
local CartDecomp     = require "Lib.CartDecomp"
local Grid           = require "Grid"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags    = require "App.Diagnostics.VlasovDiagnostics"
local DiagsImplBase  = require "App.Diagnostics.DiagnosticsImplBase"

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

-- ~~~~ Source integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
local bcRecycleDiagImpl = function()
   -- Recycling coefficient.
   local _recycleCoef = Proto(DiagsImplBase)
   function _recycleCoef:fullInit(diagApp, specIn, field, owner)
      self.field = owner.recycleCoef
      self.done  = false
   end
   function _recycleCoef:getType() return "grid" end
   function _recycleCoef:advance(tm, inFlds, outFlds) end
   -- Recycling distribution function.
   local _recycleDistF = Proto(DiagsImplBase)
   function _recycleDistF:fullInit(diagApp, specIn, field, owner)
      self.field = owner.recycleDistF
      self.done  = false
   end
   function _recycleDistF:getType() return "grid" end
   function _recycleDistF:advance(tm, inFlds, outFlds) end
   -- Recycling flux.
   local _recycleTestFlux = Proto(DiagsImplBase)
   function _recycleTestFlux:fullInit(diagApp, specIn, field, owner)
      self.field = owner.recycleTestFlux
      self.done  = false
   end
   function _recycleTestFlux:getType() return "grid" end
   function _recycleTestFlux:advance(tm, inFlds, outFlds) end
   return {recycleCoef  = _recycleCoef,
           recycleDistF = _recycleDistF,
	   recycleTestFlux = _recycleTestFlux}
end

-- .................... END OF DIAGNOSTICS ...................... --

local NeutralRecyclingBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function NeutralRecyclingBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function NeutralRecyclingBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.recycleTemp  = assert(tbl.recycleTemp, "NeutralRecyclingBC: must specify temperature using 'recycleTemp'.")
   self.recycleFrac  = assert(tbl.recycleFrac, "NeutralRecyclingBC: must specify recycling fraction using 'recycleFrac'.")
   self.recycleIonNm = assert(tbl.recycleIon,  "NeutralRecyclingBC: must specify ion species using 'recycleIon'.")

   if tbl.recycleTime then
      self.recycleTime = tbl.recycleTime
      self.scaledRecycleFrac = function(tm)
         return self.recycleFrac * (0.5*(1. + math.tanh(tm/self.recycleTime-1)))
      end
   else
      self.scaledRecycleFrac = function(tm) return self.recycleFrac end
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

function NeutralRecyclingBC:setName(nm) self.name = self.speciesName.."_"..nm end

function NeutralRecyclingBC:bcNeutralRecycling(dir, tm, idxIn, fIn, fOut)
   -- Note that bcRecycle only valid in dir parallel to B.
   -- This is checked when bc is created.
   local zIdx = idxIn
   zIdx[dir] = 1
   local numBasis = self.basis:numBasis()
   local f        = self.recycleDistF
   local rIdxr    = f:genIndexer()
   local rFPtr    = self.recycleDistF:get(1)
   f:fill(rIdxr(zIdx), rFPtr)
   for i = 1, numBasis do
      fOut[i] = 0
      fOut[i] = rFPtr[i]
   end

   self.basis:flipSign(dir, fOut:data(), fOut:data())
   self.basis:flipSign(dir+self.cdim, fOut:data(), fOut:data())
end

function NeutralRecyclingBC:initCrossSpeciesCoupling(species)
   -- Need to set saveFlux=true in the recycling ion's BC at the same boundary as this.
   local recIon = species[self.recycleIonNm]
   self.recIonBC = recIon.nonPeriodicBCs[string.gsub(self.name,self.speciesName.."_","")]
   self.recIonBC:setSaveFlux(true)
end

function NeutralRecyclingBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   self.mass = mySpecies.mass

   assert(self.bcDir==self.cdim, "NeutralRecyclingBC: recycling BC can only be used along the last/parallel configuration space dimension.")
   local bcFunc   = function(...) return self:bcNeutralRecycling(...) end
   local skinType = "flip"

   local vdir = self.bcDir+self.cdim
   self.bcSolver = Updater.Bc{
      onGrid = self.grid,   edge               = self.bcEdge,  
      cdim   = self.cdim,   skinLoop           = skinType,
      dir    = self.bcDir,  boundaryConditions = {bcFunc},   
      vdir   = vdir,      
   }

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid()

   -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
   self:createConfBoundaryGrid()

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
   -- Need to define methods to allocate fields defined on boundary grid (e.g. for diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid,self.basis:numBasis(),
                                 {distf:lowerGhost(),distf:upperGhost()},distf:getMetaData())
   end
   self.allocMoment = function(self)
      return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                 {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
   end
   self.allocVectorMoment = function(self, dim)
      return self:allocCartField(self.confBoundaryGrid, dim*self.basis:numBasis(),
                                 {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
   end

   -- For fMaxwell projection with density = 1.
   self.recycleFMaxwell = allocDistf()
   -- For projection of flux on ghosts.
   self.recycleFhat, self.scaledFhat = allocDistf(), allocDistf()
   -- For scaled projection of flux, passed to bc func.
   self.recycleDistF = allocDistf()
   -- 0th moment of fhat.
   self.recycleFhatM0 = self:allocMoment()
   -- Scaling factor for recycle distf.
   self.recycleCoef     = self:allocMoment()
   self.recycleTestFlux = self:allocMoment()

   local edgeval = self.bcEdge=="lower" and 1 or -1
   local mom
   if self.bcDir==1 then
      mom = edgeval == 1 and "M0Nvx" or "M0Pvx"
   elseif self.bcDir==3 then
      mom = edgeval == 1 and "M0Nvz" or "M0Pvz"
   end

   -- DistFuncMomentCalc updater for fMaxwell.
   self.calcFhatM0 = Updater.DistFuncMomentCalc {
      onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
      phaseBasis = self.basis,         moment     = mom,
   }

   self.recycleConfWeakDivide = Updater.CartFieldBinOp {
      onGrid    = self.confBoundaryGrid,  operation = "Divide",
      weakBasis = self.confBasis,         onGhosts  = false,
   }
   self.recycleConfPhaseWeakMultiply = Updater.CartFieldBinOp {
      onGrid     = self.boundaryGrid,  operation = "Multiply",
      weakBasis  = self.basis,         onGhosts  = false,
      fieldBasis = self.confBasis, 
   }

   local recycleSource = function (t, xn)
      local cdim, vdim = self.cdim, self.vdim
      local vt2        = self.recycleTemp/self.mass
      local v2         = 0.0
      for d = cdim+1, cdim+vdim do v2 = v2 + (xn[d])^2 end
      return 1.0 / math.sqrt(2*math.pi*vt2)^vdim * math.exp(-v2/(2*vt2))
   end

   local projectRecycleFMaxwell = Updater.ProjectOnBasis {
      onGrid = self.boundaryGrid,  evaluate = recycleSource,
      basis  = self.basis,         onGhosts = false,
   }
   local projectFluxFunc = Updater.ProjectFluxFunc {
      onGrid     = self.boundaryGrid,  edgeValue  = edgeval,
      phaseBasis = self.basis,         direction  = self.bcDir,
      confBasis  = self.confBasis,     onGhosts   = false,
   }

   projectRecycleFMaxwell:advance(0., {}, {self.recycleFMaxwell})
   projectFluxFunc:advance(0., {self.recycleFMaxwell}, {self.recycleFhat})
   self.calcFhatM0:advance(0., {self.recycleFhat}, {self.recycleFhatM0})
   self.recycleFhatM0:scale(-1.0*edgeval)

   -- Write out distf and flux.
   self.recycleFhat:write(string.format("%s_%s_%d.bp", self.name, 'recycleFhat', mySpecies.diagIoFrame),
      0., mySpecies.diagIoFrame, false)
   self.recycleFhatM0:write(string.format("%s_%s_%d.bp", self.name, 'recycleFhatM0', mySpecies.diagIoFrame),
      0., mySpecies.diagIoFrame, false)

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.boundaryPtr        = {}
      self.distfInIdxr        = distf:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = allocDistf()
         self.boundaryPtr[i]        = self.boundaryFluxFields[i]:get(1)
      end
      self.boundaryFluxRate      = allocDistf()
      self.boundaryFluxFieldPrev = allocDistf()
      self.boundaryIdxr          = self.boundaryFluxFields[1]:genIndexer()

      self.idxOut = Lin.IntVec(self.grid:ndim())

      -- Create the range needed to loop over ghosts.
      local global, globalExt, localExtRange = distf:globalRange(), distf:globalExtRange(), distf:localExtRange()
      self.ghostRng = localExtRange:intersect(self:getGhostRange(global, globalExt))
      -- Decompose ghost region into threads.
      self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRng, numSplit=self.grid:numSharedProcs()}
      self.tId              = self.grid:subGridSharedId() -- Local thread ID.

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField    = self:allocMoment()
      self.confBoundaryFieldPtr = self.confBoundaryField:get(1)
      self.confBoundaryIdxr     = self.confBoundaryField:genIndexer()
      local confGlobal        = numDensity:globalRange()
      local confGlobalExt     = numDensity:globalExtRange()
      local confLocalExtRange = numDensity:localExtRange()
      self.confGhostRange = confLocalExtRange:intersect(self:getGhostRange(confGlobal, confGlobalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self.confGhostRangeDecomp = LinearDecomp.LinearDecompRange {range=self.confGhostRange, numSplit=self.grid:numSharedProcs()}

      local jacobGeo = externalField.geo and externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                             {jacobGeo:lowerGhost(),jacobGeo:upperGhost()}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo and externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                                {jacobGeoInv:lowerGhost(),jacobGeoInv:upperGhost()}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
      end

      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         local ptrOut = qOut:get(1)
         for idx in self.ghostRangeDecomp:rowMajorIter(self.tId) do
            idx:copyInto(self.idxOut)
            qOut:fill(self.distfInIdxr(idx), ptrOut)

            -- Before operating on ghosts, store ghost values for later flux diagnostics
            self.idxOut[self.bcDir] = 1
            self.boundaryFluxFields[rkIdx]:fill(self.boundaryIdxr(self.idxOut), self.boundaryPtr[rkIdx])
            for c = 1, qOut:numComponents() do self.boundaryPtr[rkIdx][c] = ptrOut[c] end
         end
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
            onGrid    = self.confBoundaryGrid,  operation = "Multiply",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "Divide",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         self.confWeakDotProduct = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "DotProduct",
            weakBasis = self.confBasis,         onGhosts  = true,
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
         self.numDensityCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M0",
         }
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

function NeutralRecyclingBC:createCouplingSolver(species, field, externalField)
   -- Fetch and organize information about other species that recycling BCs depends on.
   local recIon = species[self.recycleIonNm]

   self.recIonBC = recIon.nonPeriodicBCs[string.gsub(self.name,self.speciesName.."_","")]
   self.bcIonM0fluxField = self.recIonBC:allocMoment() 
end

function NeutralRecyclingBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function NeutralRecyclingBC:calcCouplingMoments(tCurr, rkIdx, species)
   self.bcIonM0fluxField:scale(self.scaledRecycleFrac(tCurr))

   -- Weak divide.
   self.recycleConfWeakDivide:advance(tCurr, {self.recycleFhatM0, self.bcIonM0fluxField}, {self.recycleCoef})

   -- Weak multiply.
   self.recycleConfPhaseWeakMultiply:advance(tCurr, {self.recycleCoef, self.recycleFMaxwell}, {self.recycleDistF})

   -- Diagnostics to check flux.
   self.recycleConfPhaseWeakMultiply:advance(tCurr, {self.recycleCoef, self.recycleFhat}, {self.scaledFhat})
   self.calcFhatM0:advance(tCurr, {self.scaledFhat}, {self.recycleTestFlux})
   -- Maybe calculated integrated M0 flux here??

   self.recycleTestFlux:scale(1.0/self.recycleFrac) -- This can be written out from KineticSpecies, if necessary.
end

function NeutralRecyclingBC:advanceCrossSpeciesCoupling(tCurr, species, outIdx)
   -- Compute the 0th moment of the ion boundary flux.
   self.recIonBC.numDensityCalc:advance(tCurr, {self.recIonBC:getBoundaryFluxFields()[outIdx]}, {self.bcIonM0fluxField})
end

function NeutralRecyclingBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function NeutralRecyclingBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function NeutralRecyclingBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function NeutralRecyclingBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      local vlasovDiags  = VlasovDiags()
      local recycleDiags = bcRecycleDiagImpl()
      for nm, v in pairs(vlasovDiags) do recycleDiags[nm] = v end
      self.diagnostics = DiagsApp{implementation = recycleDiags}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

-- These are needed to recycle the VlasovDiagnostics with NeutralRecyclingBC.
function NeutralRecyclingBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                      self.boundaryFluxRate, self.boundaryFluxRate} end
function NeutralRecyclingBC:getFlucF() return self.boundaryFluxRate end

function NeutralRecyclingBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   self.bcSolver:advance(tCurr, {}, {fIn})
end

function NeutralRecyclingBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return NeutralRecyclingBC
