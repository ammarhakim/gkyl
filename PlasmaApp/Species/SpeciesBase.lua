local Proto = require "Lib.Proto"

-- empty shell species base class
local SpeciesBase = Proto()

-- functions that must be defined by subclasses
function Species:init() end
function Species:fullInit() end
function Species:setName() end
function Species:setIoMethod() end
function Species:createGrid() end
function Species:setConfBasis() end
function Species:createBasis() end
function Species:setConfGrid() end
function Species:alloc() end
function Species:setCfl() end
function Species:ndim() end
function Species:createSolver() end
function Species:createDiagnostics() end
function Species:rkStepperFields() end
function Species:initDist() end
function Species:calcCouplingMoments() end
function Species:write() end
function Species:forwardEuler() end
function Species:applyBc() end
function Species:totalSolverTime() end
function Species:momCalcTime() end
function Species:intMomCalcTime() end
function Species:totalBcTime() end
function Species:getCharge() end
function Species:getMass() end

return SpeciesBase

