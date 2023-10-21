-- Gkyl ------------------------------------------------------------------------
--
-- Species base object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- Empty shell species base class.
local SpeciesBase = Proto()

-- Functions that must be defined by subclasses.
function SpeciesBase:init(tbl) end
function SpeciesBase:fullInit(appTbl) end
function SpeciesBase:setName(nm) self.name = nm end -- Needs to be called before fullInit().
function SpeciesBase:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs, cMap) end
function SpeciesBase:setConfBasis(basis) end
function SpeciesBase:createBasis() end
function SpeciesBase:setConfGrid(grid, myadios) end
function SpeciesBase:alloc(nRkDup) end
function SpeciesBase:setCfl(cfl) end
function SpeciesBase:getNdim() return self.ndim end
function SpeciesBase:getVdim() return self.vdim end
function SpeciesBase:createSolver() end
function SpeciesBase:createCouplingSolver() end
function SpeciesBase:createDiagnostics() end
function SpeciesBase:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end
function SpeciesBase:rkStepperFields() return { nil } end
function SpeciesBase:suggestDt() end
function SpeciesBase:setDtGlobal(dtGlobal) end
function SpeciesBase:clearCFL() end
function SpeciesBase:clearMomentFlags(species) end
function SpeciesBase:initDist() end
function SpeciesBase:initCrossSpeciesCoupling(population) end
function SpeciesBase:calcCouplingMoments() end
function SpeciesBase:calcCrossCouplingMoments(tCurr, rkIdx, population) end
function SpeciesBase:write(tm) end
function SpeciesBase:writeRestart(tm) end
function SpeciesBase:readRestart(field, externalField) return 0.0 end
function SpeciesBase:advance(tCurr, population, emIn, inIdx, outIdx) return true, GKYL_MAX_DOUBLE end
function SpeciesBase:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx) end
function SpeciesBase:updateInDirection(dir, tCurr, dt, fIn, fOut) return true, GKYL_MAX_DOUBLE end
function SpeciesBase:applyBcIdx(tCurr, field, externalField, inIdx, outIdx) end
function SpeciesBase:applyBc(tCurr, field, externalField, inIdx, outIdx) end
function SpeciesBase:applyBcInitial(tCurr, field, externalField, inIdx, outIdx)
   -- This function is a temporary hack to make initial :applyBc call in PlasmaOnCartGrid work.
   self:applyBc(tCurr, field, externalField, inIdx, outIdx)
end
function SpeciesBase:getCharge() return self.charge end
function SpeciesBase:getMass() return self.mass end
function SpeciesBase:copyRk(outIdx, aIdx) end
function SpeciesBase:combineRk() end
function SpeciesBase:isEvolving() return self.evolve end
function SpeciesBase:clearTimers() end
function SpeciesBase:getTimer(timerNm)
   if self.timers then
      if self.timers[timerNm] == nil then return 0. end
      return self.timers[timerNm]
   end
   return 0.
end

return SpeciesBase

