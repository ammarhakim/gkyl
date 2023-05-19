-- Gkyl ------------------------------------------------------------------------
--
-- Maxwellian Boundary Condition for a Gyrokinetic species
-- Applied with Updater/MaxwellGhostBc
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBaseGK    = require "App.BCs.BCsBaseGK"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local Range      = require "Lib.Range"
local Lin        = require "Lib.Linalg"
local Grid       = require "Grid"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

local MaxwellianBC = Proto(BCsBaseGK)

-- Store table passed to it and defer construction to :fullInit().
function MaxwellianBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function MaxwellianBC:fullInit(mySpecies)
   self.maxwellianKind = self.tbl.maxwellianKind
   self.densityGhost   = self.tbl.densityGhost
   self.uParGhost      = self.tbl.uParGhost
   self.tempGhost      = self.tbl.tempGhost
   self.initFunc       = self.tbl.initFunc
   self.fromFile       = self.tbl.fromFile
   self.needsBoundaryTools = true
end

function MaxwellianBC:setName(nm) self.name = self.speciesName.."_"..nm end

function MaxwellianBC:createSolver(mySpecies, field, externalField)
   self.mass = mySpecies.mass
   self.charge= mySpecies.charge
   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   self:createBoundaryTools(mySpecies,field,externalField)
   self.ghostFld = self.allocDistf()
   self.phaseFieldIo = AdiosCartFieldIo {
      elemType   = mySpecies:getDistF():elemType(),
      method     = "MPI",
      writeGhost = false,
      metaData   = {polyOrder = self.basis:polyOrder(),  basisType = self.basis:id(),
                    charge    = self.charge,             mass      = self.mass,},
   }
   if not self.fromFile then
      if self.maxwellianKind == 'local' then
         self.projMaxwell = Updater.MaxwellianOnBasis{
            onGrid     = self.boundaryGrid,     confBasis = self.confBasis,
            phaseBasis = self.basis ,           mass      = self.mass,
            confGrid   = self.confBoundaryGrid, onGhosts  = false,
            usePrimMoms=true
         }
         local confProject = Updater.ProjectOnBasis {
            onGrid   = self.confBoundaryGrid,         basis    = self.confBasis,
            evaluate = function(t, xn) return 0. end, onGhosts = false
         }
         local confVec2Project = Updater.ProjectOnBasis {
            onGrid = self.confBoundaryGrid,   evaluate = function(t, xn) return 0., 0. end,
            basis  = self.confBasis,          onGhosts = false
         }
         local numDens = self:allocMoment()
         local primMoms = self:allocVectorMoment(2)
         confProject:setFunc(function(t, xn) return 1. end)
         confProject:advance(time, {}, {numDens})
         confVec2Project:setFunc(function(t,xn) return self.uParGhost(t,xn), self.tempGhost(t,xn)/self.mass end)
         confVec2Project:advance(time, {}, {primMoms})
         self.projMaxwell:advance(time,{numDens,primMoms,self.bmag,self.bmag},{self.ghostFld})
         local M0e, M0 = self:allocMoment(), self:allocMoment()
         local M0mod   = self:allocMoment()
         self.numDensityCalc:advance(0.0, {self.ghostFld}, {M0})
         confProject:setFunc(self.densityGhost)
         confProject:advance(0.0, {}, {M0e})
         self.confWeakDivide:advance(0.0, {M0, M0e}, {M0mod})
         self.phaseWeakMultiply:advance(0.0, {M0mod,self.ghostFld}, {self.ghostFld})
      elseif self.maxwellianKind == 'canonical' then
         if mySpecies.jacobPhaseFunc and self.vdim > 1 then
            local initFuncWithoutJacobian = self.initFunc
            self.initFunc = function (t, xn)
               local xconf = {}
               for d = 1, self.cdim do xconf[d] = xn[d] end
               local J = mySpecies.jacobPhaseFunc(t,xconf)
               local f = initFuncWithoutJacobian(t,xn)
               return J*f
            end
         end
         local project = Updater.ProjectOnBasis {
            onGrid   = self.boundaryGrid,
            basis    = self.basis,
            evaluate = self.initFunc,
            onGhosts = false
         }
         project:advance(time, {},{self.ghostFld})
      end
   end
end

function MaxwellianBC:createCouplingSolver(species,field,externalField)
   if self.bcKind == 'maxwellianGhost' then
      if self.fromFile then
         local tm, fr = self.phaseFieldIo:read(self.ghostFld, self.fromFile)
      else
         if self.maxwellianKind == 'canonical' then
            local ionName = nil
            for nm, s in lume.orderedIter(species) do
               if 0.0 < s.charge then ionName = nm end
            end

            local numDens = self:allocMoment()
            if self.species.charge < 0.0 then
               local bcName = self.name:gsub(self.speciesName .. "_","")
               local numDensScaleTo = species[ionName]['nonPeriodicBCs'][bcName]:allocMoment()
               self.numDensityCalc:advance(0.0, {self.ghostFld}, {numDens})
               species[ionName]['nonPeriodicBCs'][bcName].numDensityCalc:advance(0.0,
                  {species[ionName]['nonPeriodicBCs'][bcName].ghostFld}, {numDensScaleTo})

               local M0mod = self:allocMoment()
               self.confWeakDivide:advance(0.0,{numDens, numDensScaleTo} , {M0mod})
               self.phaseWeakMultiply:advance(0.0, {M0mod,self.ghostFld}, {self.ghostFld})
            end
            self.numDensityCalc:advance(0.0, {self.ghostFld}, {numDens})
            numDens:write(string.format("%s_M0_%d.bp", self.name, 0), 0.0, 0, false)
         end
         -- Multiply by conf space Jacobian Now
         self.phaseWeakMultiply:advance(0, {self.ghostFld, self.jacobGeo}, {self.ghostFld})
         self.phaseFieldIo:write(self.ghostFld, string.format("%s_Maxwellian_%d.bp", self.name, 0), 0., 0, false)
      end
   end
end

function MaxwellianBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function MaxwellianBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function MaxwellianBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function MaxwellianBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function MaxwellianBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with MaxwellianBC.
function MaxwellianBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function MaxwellianBC:getFlucF() return self.boundaryFluxRate end

function MaxwellianBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   fIn:copyRangeToRange(self.ghostFld, self.myGlobalGhostRange,self.ghostFld:localRange())
end

function MaxwellianBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return MaxwellianBC
