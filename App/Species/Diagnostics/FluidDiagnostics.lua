-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics specific to FluidSpecies.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto         = require "Lib.Proto"
local DiagsImplBase = require "App.Species.Diagnostics.DiagnosticsImplBase"
local DataStruct    = require "DataStruct"
local Updater       = require "Updater"

-- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
local FluidDiag_MomSq = Proto(DiagsImplBase)
function FluidDiag_MomSq:fullInit(diagApp, specIn)
   self.field    = specIn:allocVectorMoment(specIn.nMoments)
   self.updaters = specIn.weakMultiply
   self.done     = false
end
function FluidDiag_MomSq:getDependencies() return {} end
function FluidDiag_MomSq:getType() return "grid" end
function FluidDiag_MomSq:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   local momIn  = specIn:getMoments(1)
   self.updaters:advance(tm, {momIn, momIn}, {self.field})
end

-- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
local FluidDiag_intMom = Proto(DiagsImplBase)
function FluidDiag_intMom:fullInit(diagApp, specIn)
   self.field    = DataStruct.DynVector { numComponents = specIn.nMoments }
   self.updaters = specIn.volIntegral.compsN
   self.done     = false
end
function FluidDiag_intMom:getDependencies() return {} end
function FluidDiag_intMom:getType() return "integrated" end
function FluidDiag_intMom:advance(tm, inFlds, outFlds)
   local specIn = inFlds[1]
   self.updaters:advance(tm, {specIn:rkStepperFields()[1]}, {self.field})
end

return {
  MomSq  = FluidDiag_MomSq,
  intMom = FluidDiag_intMom
}
