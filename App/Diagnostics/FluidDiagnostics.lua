-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics specific to FluidSpecies.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto         = require "Lib.Proto"
local DiagsImplBase = require "App.Diagnostics.DiagnosticsImplBase"
local DataStruct    = require "DataStruct"
local Updater       = require "Updater"

local implementation = function() 
   -- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
   local _MomSq = Proto(DiagsImplBase)
   function _MomSq:fullInit(diagApp, specIn)
      self.field    = specIn:allocVectorMoment(specIn.nMoments)
      self.done     = false
   end
   function _MomSq:getType() return "grid" end
   function _MomSq:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local momIn  = specIn:getMoments(1)
      specIn.weakMultiply:advance(tm, {momIn, momIn}, {self.field})
   end
   
   -- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
   local _intMom = Proto(DiagsImplBase)
   function _intMom:fullInit(diagApp, specIn)
      self.field    = DataStruct.DynVector { numComponents = specIn.nMoments }
      self.updaters = specIn.volIntegral.compsN
      self.done     = false
   end
   function _intMom:getType() return "integrated" end
   function _intMom:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      self.updaters:advance(tm, {specIn:rkStepperFields()[1]}, {self.field})
   end

   return {
     MomSq  = _MomSq,
     intMom = _intMom
   }
end

return implementation
