-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic sheath BC.
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

local GkSheathBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GkSheathBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkSheathBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.phiWallIn = tbl.phiWall
   if self.phiWallIn then assert(type(self.phiWallIn)=="table" or type(self.phiWallIn)=="function", "GkSheathBC: phiWall must be a function or a table with 'spatialDependence', and optionally 'timeDependence', functions.") end

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

function GkSheathBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GkSheathBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   -- Create the boundary grid and other boundary tools.
   BCtools.createBoundaryTools(mySpecies, field, externalField, self)

   local charge, mass = mySpecies.charge, mySpecies.mass

   -- Create field and function for calculating wall potential according to user-provided function.
   self.phiWallFld = self:allocCartField(field.grid, field.basis:numBasis(), {1,1},
                                         {polyOrder=field.basis:polyOrder(), basisType=field.basis:id()},false)
   self.phiWallFld:clear(0.)
   if self.phiWallIn then
      local projPhiWall = Updater.EvalOnNodes {
         onGrid = field.grid,   evaluate = function(t,xn) return 1. end,
         basis  = field.basis,  onGhosts = true,
      }
      if type(self.phiWallIn) == 'function' then
         self.setPhiWall = projPhiWall
         self.setPhiWall:setFunc(self.phiWallIn)
         self.setPhiWall:advance(0.0, {}, {self.phiWallFld})
      elseif type(self.phiWallIn) == 'table' then
         projPhiWall:setFunc(function(t,xn) return self.phiWallIn["spatialDependence"](xn) end)
         if self.phiWallIn["timeDependence"] then
            self.phiWallFldInit = self:allocCartField(field.grid, field.basis:numBasis(), {1,1},
                                                      {polyOrder=field.basis:polyOrder(), basisType=field.basis:id()},false)
            projPhiWall:advance(self.phiWallIn["timeDependence"](0.), {}, {self.phiWallFldInit})
            self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds)
               self.phiWallFld:combine(self.phiWallIn["timeDependence"](tCurr), self.phiWallFldInit)
            end}
         else
            -- No time dependence
            projPhiWall:advance(0., {}, {self.phiWallFld})
            self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
         end
      end
   else
      self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
   end
   self.phiWallFld:sync(false)

   self.bcSolver = Updater.GkSheathBc{
      onGrid  = self.grid,   edge    = self.bcEdge,  
      cdim    = self.cdim,   basis   = self.basis,
      dir     = self.bcDir,  phiWall = self.phiWallFld,
      charge  = charge,      mass    = mass,
      onField = mySpecies:rkStepperFields()[1],
   }

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

function GkSheathBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function GkSheathBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GkSheathBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GkSheathBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GkSheathBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with GkSheathBC.
function GkSheathBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function GkSheathBC:getFlucF() return self.boundaryFluxRate end

function GkSheathBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   self.setPhiWall:advance(tCurr, {}, {self.phiWallFld}) -- Compute wall potential if needed.

   local phi = field:rkStepperFields()[inIdx].phi         -- Get the current plasma potential.
   local fIn = mySpecies:rkStepperFields()[outIdx] 

   self.bcSolver:advance(tCurr, {phi}, {fIn})
end

function GkSheathBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return GkSheathBC
