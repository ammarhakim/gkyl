-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a Gyrokinetic species, i.e. those that can be
-- applied with Updater/Bc.lua using just a function (e.g. bcAbsorb)
-- and no additional setup.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Proto      = require "Lib.Proto"
local BCtools    = require "App.BCs.GkBCtools"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"

local GkBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GkBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkBasicBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind = assert(tbl.kind, "GkBasicBC: must specify the type of BC in 'kind'.")

   self.evolve = xsys.pickBool(tbl.evolve, false) 

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function GkBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GkBasicBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   -- Create the boundary grid and other boundary tools.
   BCtools.createBoundaryTools(mySpecies, field, externalField, self)

   self.bcBuffer = self:allocDistf() -- Buffer used by BasicBc updater.

   local bcFunc, skinType
   if self.bcKind == "copy" or self.bcKind == "absorb" or self.bcKind == "reflect" then
      self.bcSolver = Updater.BasicBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   basis  = self.basis,
         dir     = self.bcDir,  bcType = self.bcKind,
         onField = mySpecies:rkStepperFields()[1],
      }
      self.bcSolverAdvance = function(tm, inFlds, outFlds)
         self.bcSolver:advance(tm, {inFlds[1]}, outFlds)
      end
   else
      assert(false, "GkBasicBC: bcKind not supported.")
   end

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
end

function GkBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function GkBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GkBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GkBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GkBasicBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with GkBasicBC.
function GkBasicBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function GkBasicBC:getFlucF() return self.boundaryFluxRate end

function GkBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   self.bcSolverAdvance(tCurr, {self.bcBuffer}, {fIn})
end

function GkBasicBC:getBoundaryFluxFields() return self.boundaryFluxFields end

-- ................... Classes meant as aliases to simplify input files ...................... --
local GkAbsorbBC = Proto(GkBasicBC)
function GkAbsorbBC:fullInit(mySpecies)
   self.tbl.kind  = "absorb"
   GkAbsorbBC.super.fullInit(self, mySpecies)
end

local GkReflectBC = Proto(GkBasicBC)
function GkReflectBC:fullInit(mySpecies)
   self.tbl.kind  = "reflect"
   GkReflectBC.super.fullInit(self, mySpecies)
end

local GkCopyBC = Proto(GkBasicBC)
function GkCopyBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   GkCopyBC.super.fullInit(self, mySpecies)
end

local GkZeroFluxBC = Proto()
function GkZeroFluxBC:init(tbl)
   self.tbl      = tbl
   self.tbl.kind = "zeroFlux"
end

-- ................... End of GkBasicBC alias classes .................... --

return {GkBasic    = GkBasicBC,
        GkAbsorb   = GkAbsorbBC,
        GkCopy     = GkCopyBC,
        GkReflect  = GkReflectBC,
        GkZeroFlux = GkZeroFluxBC,}
