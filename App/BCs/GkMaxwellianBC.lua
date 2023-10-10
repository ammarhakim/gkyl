-- Gkyl ------------------------------------------------------------------------
--
-- Maxwellian Boundary Condition for a Gyrokinetic species
-- Applied with Updater/MaxwellGhostBc
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Proto      = require "Lib.Proto"
local lume       = require "Lib.lume"
local BCtools    = require "App.BCs.GkBCtools"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

local GkMaxwellianBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GkMaxwellianBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkMaxwellianBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.densityGhost   = tbl.density
   self.uParGhost      = tbl.driftSpeed or function(t,xn) return 0. end
   self.tempGhost      = tbl.temperature
   self.maxwellianKind = tbl.kind or 'local'
   self.profileIn      = tbl.profile
   self.fromFile       = tbl.fromFile

   assert((self.densityGhost and self.tempGhost) or self.fromFile or self.profileIn, "GkMaxwellianBC: must specify the density and temperature of the Maxwellian, or the 'fromFile' option to read the Maxwellian from a file, or the phase-space profile with 'profile'.")
   assert(self.maxwellianKind=='local' or self.maxwellianKind=='canonical', "GkMaxwellianBC: 'kind' must be local or canonical.")

   self.mass, self.charge = mySpecies.mass, mySpecies.charge

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function GkMaxwellianBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GkMaxwellianBC:createSolver(mySpecies, field, externalField)
   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   -- Create the boundary grid and other boundary tools.
   BCtools.createBoundaryTools(mySpecies, field, externalField, self)

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then

      -- Create boundary tools for saving fluxes.
      BCtools.createFluxTools(mySpecies, field, externalField, self)

      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         -- Create boundary tools for diagnostics.
         BCtools.createDiagnosticTools(mySpecies, field, externalField, self)
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end

   self.ghostFld = self:allocDistf()
   -- The following range is a bit of a hack to get the copyRangeToRange in :advance
   -- to only take place in the MPI process that owns the ghost range.
   self.ghostFldLocalGlobalGhostRange = self.myGlobalGhostRange:volume() < 1 and self.myGlobalGhostRange or self.ghostFld:localRange()

   if self.fromFile then
      local phaseFieldIo = AdiosCartFieldIo {
         elemType   = mySpecies:getDistF():elemType(),
         method     = "MPI",  writeGhost = false,
         metaData   = {polyOrder = mySpecies.basis:polyOrder(),  basisType = mySpecies.basis:id(),
                       charge    = mySpecies.charge,             mass      = mySpecies.mass,},
      }

      local tm, fr = phaseFieldIo:read(self.ghostFld, self.fromFile)
   else
      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      local confBoundaryField = self.confBoundaryField or self:allocMoment()
      local bmag = externalField.geo.bmag
      self.bmag = self.bmag or self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmag:getMetaData())
      confBoundaryField = self:evalOnConfBoundary(bmag, confBoundaryField)
      self.bmag:copy(self:evalOnConfBoundary(bmag, confBoundaryField))

      self.numDensityCalc = self.numDensityCalc or Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
         phaseBasis = self.basis,         gkfacs    = {self.mass, self.bmag},
         moment     = "GkM0", -- GkM0 = < f >
      }

      if self.maxwellianKind == 'local' then
         local projMaxwell = Updater.MaxwellianOnBasis{
            onGrid      = self.boundaryGrid,      confBasis = self.confBasis,
            phaseBasis  = self.basis ,            mass      = self.mass,
            confGrid    = self.confBoundaryGrid,  onGhosts  = false,
            usePrimMoms = true
         }
         local confProject = Updater.ProjectOnBasis {
            onGrid = self.confBoundaryGrid,  evaluate = function(t, xn) return 0. end,
            basis  = self.confBasis,         onGhosts = false
         }
         local confVec2Project = Updater.ProjectOnBasis {
            onGrid = self.confBoundaryGrid,  evaluate = function(t, xn) return 0., 0. end,
            basis  = self.confBasis,         onGhosts = false
         }
         local numDens, primMoms = self:allocMoment(), self:allocVectorMoment(2)
   
         local jacobGeo = externalField.geo.jacobGeo
         if jacobGeo then
            self.jacobGeo = self.jacobGeo or self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
            self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo, confBoundaryField))
         end
   
         confProject:setFunc(function(t, xn) return 1. end)
         confProject:advance(time, {}, {numDens})
         confVec2Project:setFunc(function(t,xn) return self.uParGhost(t,xn), self.tempGhost(t,xn)/self.mass end)
         confVec2Project:advance(time, {}, {primMoms})
         projMaxwell:advance(time,{numDens,primMoms,self.bmag,self.bmag},{self.ghostFld})
   
         local M0e, M0, M0mod = self:allocMoment(), self:allocMoment(), self:allocMoment()
         local confWeakDivide = self.confWeakDivide or Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "Divide",
            onRange   = confBoundaryField:localRange(),  onGhosts = false,
         }
         local phaseWeakMultiply = self.phaseWeakMultiply or Updater.CartFieldBinOp {
            onGrid    = self.boundaryGrid,  operation  = "Multiply",
            weakBasis = self.basis,         fieldBasis = self.confBasis,
            onGhosts  = true,
         }
         self.numDensityCalc:advance(0.0, {self.ghostFld}, {M0})
         confProject:setFunc(self.densityGhost)
         confProject:advance(0.0, {}, {M0e})
         confWeakDivide:advance(0.0, {M0, M0e}, {M0mod})
         phaseWeakMultiply:advance(0.0, {M0mod,self.ghostFld}, {self.ghostFld})

      else
         -- This option is intended for Canonical Maxwellians, which are currently defined via a
         -- function in the input file, although in principle one could also specify other types
         -- of profiles.
         local initFunc = self.profileIn
         if mySpecies.jacobPhaseFunc and self.vdim > 1 then
            initFunc = function(t, xn)
               local xconf = {}
               for d = 1, self.cdim do xconf[d] = xn[d] end
               local J = mySpecies.jacobPhaseFunc(t,xconf)
               local f = self.profileIn(t,xn)
               return J*f
            end
         end
         local project = Updater.ProjectOnBasis {
            onGrid = self.boundaryGrid,  evaluate = initFunc,
            basis  = self.basis,         onGhosts = false
         }
         project:advance(time, {},{self.ghostFld})
      end
   end
end

function GkMaxwellianBC:createCouplingSolver(species,field,externalField)
   if not self.fromFile then
      local mySpecies = species[self.speciesName]
      local basis, confBasis = mySpecies.basis, mySpecies.confBasis
      local phaseFieldIo = AdiosCartFieldIo {
         elemType   = field:rkStepperFields()[1].phi:elemType(),
         method     = "MPI",  writeGhost = false,
         metaData   = {polyOrder = basis:polyOrder(),  basisType = basis:id(),
                       charge    = self.charge,             mass      = self.mass,},
      }

      local phaseWeakMultiply = self.phaseWeakMultiply or Updater.CartFieldBinOp {
         onGrid    = self.boundaryGrid,  operation  = "Multiply",
         weakBasis = basis,         fieldBasis = confBasis,
         onGhosts  = true,
      }

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      local confBoundaryField = self.confBoundaryField or self:allocMoment()
      local bmag = externalField.geo.bmag
      self.bmag = self.bmag or self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmag:getMetaData())
      confBoundaryField = self:evalOnConfBoundary(bmag, confBoundaryField)
      self.bmag:copy(self:evalOnConfBoundary(bmag, confBoundaryField))
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self.jacobGeo or self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo, confBoundaryField))
      end

      if self.maxwellianKind == 'canonical' then
         local confBoundaryField = self.confBoundaryField or self:allocMoment()
         local confWeakDivide = self.confWeakDivide or Updater.CartFieldBinOp {
            weakBasis = confBasis,  operation = "Divide",
            onRange   = confBoundaryField:localRange(),  onGhosts = false,
         }

         local ionName = nil
         for nm, s in lume.orderedIter(species) do
            if 0.0 < s.charge then ionName = nm end
         end

         local numDens = self:allocMoment()
         if self.charge < 0.0 then
            -- Scale the electron Maxwellian to have the same density as the ions.
            local bcName = self.name:gsub(self.speciesName .. "_","")
            local numDensScaleTo = species[ionName]['nonPeriodicBCs'][bcName]:allocMoment()
            self.numDensityCalc:advance(0.0, {self.ghostFld}, {numDens})
            species[ionName]['nonPeriodicBCs'][bcName].numDensityCalc:advance(0.0,
               {species[ionName]['nonPeriodicBCs'][bcName].ghostFld}, {numDensScaleTo})

            local M0mod = self:allocMoment()
            confWeakDivide:advance(0.0,{numDens, numDensScaleTo} , {M0mod})
            phaseWeakMultiply:advance(0.0, {M0mod,self.ghostFld}, {self.ghostFld})
         end
      end
      -- Multiply by conf space Jacobian.
      if self.ghostFld then  -- If parallelizeSpecies, only species owned by this MPI process has self.ghostFld.
         phaseWeakMultiply:advance(0, {self.ghostFld, self.jacobGeo}, {self.ghostFld})
         phaseFieldIo:write(self.ghostFld, string.format("%s_Maxwellian_%d.bp", self.name, 0), 0., 0, false)
      end
   end
end

function GkMaxwellianBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function GkMaxwellianBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GkMaxwellianBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GkMaxwellianBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GkMaxwellianBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GkDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

-- These are needed to recycle the GkDiagnostics with GkMaxwellianBC.
function GkMaxwellianBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                  self.boundaryFluxRate, self.boundaryFluxRate} end
function GkMaxwellianBC:getFlucF() return self.boundaryFluxRate end

function GkMaxwellianBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   fIn:copyRangeToRange(self.ghostFld, self.myGlobalGhostRange, self.ghostFldLocalGlobalGhostRange)
end

function GkMaxwellianBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return GkMaxwellianBC
