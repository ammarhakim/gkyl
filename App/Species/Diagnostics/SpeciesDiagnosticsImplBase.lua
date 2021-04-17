-- Gkyl ------------------------------------------------------------------------
--
-- Base object for the implementation of species diagnostics.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell species base class.
local SpeciesDiagnosticsImplBase = Proto()

-- Functions that must be defined by subclasses.
function SpeciesDiagnosticsImplBase:fullInit(diagApp, thisSpecies) end
function SpeciesDiagnosticsImplBase:calc(time, thisSpecies) end
function SpeciesDiagnosticsImplBase:getDependencies() return {} end
function SpeciesDiagnosticsImplBase:getType() end
function SpeciesDiagnosticsImplBase:advance(time, inFlds, outFlds) end

return SpeciesDiagnosticsImplBase
