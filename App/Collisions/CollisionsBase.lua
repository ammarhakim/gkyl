local Proto = require "Lib.Proto"

-- Empty shell collisions base class.
local CollisionsBase = Proto()

function CollisionsBase:init(tbl) self.tbl = tbl end
function CollisionsBase:setName(nm)
   self.name = self.speciesName.."_"..nm
   self.collNm = nm
end
function CollisionsBase:createSolver(mySpecies, externalField) end
function CollisionsBase:createCouplingSolver(species, field, externalField) end
function CollisionsBase:setSpeciesName(nm) self.speciesName = nm end
function CollisionsBase:fullInit(speciesTbl) end
function CollisionsBase:createDiagnostics(mySpecies, field) return nil end
function CollisionsBase:calcCouplingMoments(tCurr, rkIdx, species) end
function CollisionsBase:calcCrossCouplingMoments(tCurr, rkIdx, population) end
function CollisionsBase:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx) end
function CollisionsBase:setCfl(cfl) self.cfl = cfl end
function CollisionsBase:clearTimers() 
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
end
function CollisionsBase:getTimer(timerNm)
   if self.timers then
      if self.timers[timerNm] == nil then return 0. end
      return self.timers[timerNm]
   end
   return 0.
end

return CollisionsBase
