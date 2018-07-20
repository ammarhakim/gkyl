-- Gkyl ------------------------------------------------------------------------
--
-- Species base object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"

-- empty shell species base class
local SpeciesBase = Proto()

-- functions that must be defined by subclasses
function SpeciesBase:init(tbl) end
function SpeciesBase:fullInit(appTbl) end
function SpeciesBase:setName(nm) end
function SpeciesBase:setIoMethod(ioMethod) end
function SpeciesBase:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs, cMap) end
function SpeciesBase:setConfBasis(basis) end
function SpeciesBase:createBasis() end
function SpeciesBase:setConfGrid(grid) end
function SpeciesBase:alloc(nRkDup) end
function SpeciesBase:setCfl(cfl) end
function SpeciesBase:getNdim() end
function SpeciesBase:createSolver() end
function SpeciesBase:createDiagnostics() end
function SpeciesBase:rkStepperFields() end
function SpeciesBase:initDist() end
function SpeciesBase:calcCouplingMoments() end
function SpeciesBase:write(tm) end
function SpeciesBase:writeRestart(tm) end
function SpeciesBase:readRestart() return 0.0 end
function SpeciesBase:forwardEuler(tCurr, dt, fIn, emIn, fOut) end
function SpeciesBase:applyBc(tCurr, dt, fld) end
function SpeciesBase:totalSolverTime() end
function SpeciesBase:momCalcTime() end
function SpeciesBase:intMomCalcTime() end
function SpeciesBase:totalBcTime() end
function SpeciesBase:getCharge() end
function SpeciesBase:getMass() end
function SpeciesBase:copyRk() end
function SpeciesBase:combineRk() end

return SpeciesBase

