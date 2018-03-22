local Proto = require "Lib.Proto"

-- empty shell species base class
local SpeciesBase = Proto()

-- functions that must be defined by subclasses
function SpeciesBase:init() end
function SpeciesBase:fullInit() end
function SpeciesBase:setName() end
function SpeciesBase:setIoMethod() end
function SpeciesBase:createGrid() end
function SpeciesBase:setConfBasis() end
function SpeciesBase:createBasis() end
function SpeciesBase:setConfGrid() end
function SpeciesBase:alloc() end
function SpeciesBase:setCfl() end
function SpeciesBase:getNdim() end
function SpeciesBase:createSolver() end
function SpeciesBase:createDiagnostics() end
function SpeciesBase:rkStepperFields() end
function SpeciesBase:initDist() end
function SpeciesBase:calcCouplingMoments() end
function SpeciesBase:write() end
function SpeciesBase:forwardEuler() end
function SpeciesBase:applyBc() end
function SpeciesBase:totalSolverTime() end
function SpeciesBase:momCalcTime() end
function SpeciesBase:intMomCalcTime() end
function SpeciesBase:totalBcTime() end
function SpeciesBase:getCharge() end
function SpeciesBase:getMass() end

return SpeciesBase

