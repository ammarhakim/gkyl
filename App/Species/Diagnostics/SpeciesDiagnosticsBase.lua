-- Gkyl ------------------------------------------------------------------------
--
-- Species diagnostics base object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell species base class.
local SpeciesDiagnosticsBase = Proto()

-- Functions that must be defined by subclasses.
function SpeciesDiagnosticsBase:init(selfSpecies) end
function SpeciesDiagnosticsBase:write() end

return SpeciesDiagnosticsBase
