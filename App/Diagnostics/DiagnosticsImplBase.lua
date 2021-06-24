-- Gkyl ------------------------------------------------------------------------
--
-- Base object for the implementation of species diagnostics.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell species base class.
local DiagnosticsImplBase = Proto()

-- Functions that must be defined by subclasses.
function DiagnosticsImplBase:fullInit(diagApp, thisSpecies, field, diagAppOwner) end
function DiagnosticsImplBase:getDependencies() return {} end
function DiagnosticsImplBase:getType() end
function DiagnosticsImplBase:setGetF(perturbed, owner)
   if perturbed then
      self.getF = function() return owner:getFlucF() end
   else
      self.getF = function() return owner:rkStepperFields()[1] end
   end
end
function DiagnosticsImplBase:advance(time, inFlds, outFlds) end

return DiagnosticsImplBase
