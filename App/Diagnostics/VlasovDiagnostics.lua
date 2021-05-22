-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for VlasovSpecies.
--
-- We assume that integrated diagnostics can depend on grid diagnostics, but
-- not the other way around. We also assume that the integrated diagnostics are
-- always computed before the field diagnostics, since they may be computed more
-- frequently. This allows us to reset the state in calcIntegratedDiagostics.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto         = require "Lib.Proto"
local DiagsImplBase = require "App.Diagnostics.DiagnosticsImplBase"
local DataStruct    = require "DataStruct"
local Updater       = require "Updater"

local implementation = function()
   -- ~~~~ Number density ~~~~~~~~~~~~~~~~~~~~~~
   local _M0 = Proto(DiagsImplBase)
   function _M0:fullInit(diagApp, specIn, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.numDensityCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _M1i = Proto(DiagsImplBase)
   function _M1i:fullInit(diagApp, specIn, owner)
      self.field = owner:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _M1i:getType() return "grid" end
   function _M1i:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      local fIn    = self.owner:rkStepperFields()[1]
      specIn.momDensityCalc:advance(tm, {fIn}, {self.field})
      specIn.divideByJacobGeo(tm, self.field)
   end

   -- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
   local _intM0 = Proto(DiagsImplBase)
   function _intM0:fullInit(diagApp, specIn, owner)
      self.field    = DataStruct.DynVector { numComponents = 1 }
      self.done     = false
   end
   function _intM0:getDependencies() return {"M0"} end
   function _intM0:getType() return "integrated" end
   function _intM0:advance(tm, inFlds, outFlds)
      local specIn, diags = inFlds[1], inFlds[2]
      local M0 = diags["M0"].field
      specIn.volIntegral.comps1:advance(tm, {M0}, {self.field})
   end

   return {
      M0     = _M0,
      M1i    = _M1i,
      intM0     = _intM0,
   }
end

return implementation

