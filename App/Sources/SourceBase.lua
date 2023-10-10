-- Gkyl ------------------------------------------------------------------------
--
-- Source base object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell source base class.
local SourceBase = Proto()

-- Functions that must be defined by subclasses.
function SourceBase:init(tbl) self.tbl = tbl end
function SourceBase:fullInit(speciesTbl) end
function SourceBase:setName(nm) self.name = nm end
function SourceBase:setSpeciesName(nm) self.speciesName = nm end
function SourceBase:setConfBasis(basis) self.confBasis = basis end
function SourceBase:setConfGrid(grid) self.confGrid = grid end
function SourceBase:initCrossSpeciesCoupling(population) end
function SourceBase:createSolver(thisSpecies, extField) end
function SourceBase:createCouplingSolver(population, field, extField) end
function SourceBase:advance(tCurr, fIn, species, fRhsOut) end
function SourceBase:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end
function SourceBase:createDiagnostics(thisSpecies, momTable) end
function SourceBase:writeDiagnosticIntegratedMoments(tm, frame) end
function SourceBase:write(tm, frame) end
function SourceBase:clearTimers() 
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
end
function SourceBase:getTimer(timerNm)
   if self.timers then
      if self.timers[timerNm] == nil then return 0. end
      return self.timers[timerNm]
   end
   return 0.
end

return SourceBase
