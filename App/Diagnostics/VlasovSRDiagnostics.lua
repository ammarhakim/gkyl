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
   function _M0:fullInit(diagApp, specIn, field, owner)
      self.field   = owner:allocMoment()
      self.done = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      specIn:getNumDensity()
      self.field:clear(0.0)
      self.field:accumulate(1.0, specIn.numDensity)
   end

   -- ~~~~ Momentum density ~~~~~~~~~~~~~~~~~~~~~~
   local _M1i = Proto(DiagsImplBase)
   function _M1i:fullInit(diagApp, specIn, field, owner)
      self.field = owner:allocVectorMoment(specIn.vdim)
      self.done  = false
   end
   function _M1i:getType() return "grid" end
   function _M1i:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      specIn:getMomDensity()
      self.field:clear(0.0)
      self.field:accumulate(1.0, specIn.momDensity)
   end

   return {
      M0        = _M0,
      M1i       = _M1i,
   }
end

return implementation

