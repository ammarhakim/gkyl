-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a Vlasov species, i.e. those that can be
-- applied with Updater/Bc.lua using just a function (e.g. bcAbsorb, bcOpen)
-- and no additional setup.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBaseVlasov     = require "App.BCs.BCsBaseVlasov"
local DataStruct  = require "DataStruct"
local Updater     = require "Updater"
local Mpi         = require "Comm.Mpi"
local Proto       = require "Lib.Proto"
local Time        = require "Lib.Time"
local Range       = require "Lib.Range"
local Lin         = require "Lib.Linalg"
local Grid        = require "Grid"
local DiagsApp    = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags = require "App.Diagnostics.VlasovDiagnostics"
local xsys        = require "xsys"

local VlasovBasicBC = Proto(BCsBaseVlasov)

-- Store table passed to it and defer construction to :fullInit().
function VlasovBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VlasovBasicBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   if tbl.kind=="function" or tbl.bcFunction then
      assert(type(tbl.bcFunction)=="function", "VlasovBasicBC: bcFunction must be a function.")
      self.bcKind   = "function"
      self.bcFuncIn = assert(tbl.bcFunction, "VlasovBasicBC: must specify the BC function in 'bcFunc' when using 'function' BC kind.")
      self.feedback = xsys.pickBool(tbl.feedback, false) 
      self.evolve   = xsys.pickBool(tbl.evolve, self.feedback) 
   else
      self.bcKind = assert(tbl.kind, "VlasovBasicBC: must specify the type of BC in 'kind'.")
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

function VlasovBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function VlasovBasicBC:bcOpen(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   self.basis:flipSign(dir, fIn:data(), fOut:data())
end

function VlasovBasicBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()
   self:createBoundaryTools(mySpecies, field, externalField)

   self.bcBuffer = self.allocDistf() -- Buffer used by BasicBc updater.

   local bcFunc, skinType
   if self.bcKind == "copy" or self.bcKind == "absorb" or self.bcKind == "reflect" then
      self.bcSolver = Updater.BasicBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   basis  = self.basis,
         dir     = self.bcDir,  bcType = self.bcKind,
         onField = mySpecies:rkStepperFields()[1],
      }
   else
      -- g2, to be deleted.
      if self.bcKind == "open" then
         bcFunc   = function(...) return self:bcOpen(...) end
         skinType = "pointwise"
      elseif self.bcKind == "function" then
         bcFunc   = function(...) return self:bcCopy(...) end
         skinType = "pointwise"
      else
         assert(false, "VlasovBasicBC: BC kind not recognized.")
      end

      local vdir = self.bcDir+self.cdim

      self.bcSolver = Updater.Bc{
         onGrid   = self.grid,   edge               = self.bcEdge,  
         cdim     = self.cdim,   boundaryConditions = {bcFunc},   
         dir      = self.bcDir,  evaluate           = self.bcFuncIn,
         vdir     = vdir,        evolveFn           = self.evolve,
         skinLoop = skinType,    feedback           = self.feedback,
         basis    = self.basis,  confBasis          = self.confBasis,
         advanceArgs = {{mySpecies:rkStepperFields()[1]}, {mySpecies:rkStepperFields()[1]}},
      }
   end

end

function VlasovBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function VlasovBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function VlasovBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function VlasovBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function VlasovBasicBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the VlasovDiagnostics with VlasovBasicBC.
function VlasovBasicBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                 self.boundaryFluxRate, self.boundaryFluxRate} end
function VlasovBasicBC:getFlucF() return self.boundaryFluxRate end

function VlasovBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]

   self.bcSolver:advance(tCurr, {self.bcBuffer}, {fIn})
end

function VlasovBasicBC:getBoundaryFluxFields() return self.boundaryFluxFields end

-- ................... Classes meant as aliases to simplify input files ...................... --
local VlasovAbsorbBC = Proto(VlasovBasicBC)
function VlasovAbsorbBC:fullInit(mySpecies)
   self.tbl.kind  = "absorb"
   VlasovAbsorbBC.super.fullInit(self, mySpecies)
end

local VlasovReflectBC = Proto(VlasovBasicBC)
function VlasovReflectBC:fullInit(mySpecies)
   self.tbl.kind  = "reflect"
   VlasovReflectBC.super.fullInit(self, mySpecies)
end

local VlasovCopyBC = Proto(VlasovBasicBC)
function VlasovCopyBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   VlasovCopyBC.super.fullInit(self, mySpecies)
end

local VlasovOpenBC = Proto(VlasovBasicBC)
function VlasovOpenBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   VlasovOpenBC.super.fullInit(self, mySpecies)
end

local VlasovZeroFluxBC = Proto()
function VlasovZeroFluxBC:init(tbl)
   self.tbl      = tbl
   self.tbl.kind = "zeroFlux"
end
-- ................... End of VlasovBasicBC alias classes .................... --

return {VlasovBasic    = VlasovBasicBC,
        VlasovAbsorb   = VlasovAbsorbBC,
        VlasovCopy     = VlasovCopyBC,
        VlasovOpen     = VlasovOpenBC,
        VlasovReflect  = VlasovReflectBC,
        VlasovZeroFlux = VlasovZeroFluxBC}
